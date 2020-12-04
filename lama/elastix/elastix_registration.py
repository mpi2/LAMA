from logzero import logger as logging
from os.path import join, isdir, splitext, basename, relpath
import subprocess
import os
import shutil
from collections import defaultdict
from pathlib import Path
import yaml
import SimpleITK as sitk

from lama.elastix.folding import unfold_bsplines
from lama import common
from lama.elastix import ELX_TRANSFORM_NAME, RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR

RESOLUTION_TP_PREFIX = 'TransformParameters.0.R'
FULL_STAGE_TP_FILENAME = 'TransformParameters.0.txt'


class ElastixRegistration(object):

    def __init__(self, elxparam_file: Path,
                 movdir: Path,
                 stagedir: Path,
                 filetype: str,
                 threads: int = None,
                 fixed_mask=None,
                 ):
        """

        Parameters
        ----------
        elxparam_file
            The elastix parameter file (in elastix format)
        movdir
            Can be a directory of images or a path to a single image
        stagedir
            The output directory
        filetype
            Usually ?
        threads
            How many threads to maximally use
        fixed_mask
            The binary mask for fixed image
        """

        self.elxparam_file = elxparam_file
        self.movdir = movdir
        self.stagedir = stagedir
        self.fixed_mask = fixed_mask
        self.filetype = filetype
        self.threads = threads
        self.rename_output = True  # Bodge for pairwise reg, or we end up filling all the disks


    def make_average(self, out_path):
        """
        Create an average of the the input embryo volumes.
        This will search subfolders for all the registered volumes within them
        """
        vols = common.get_file_paths(self.stagedir, ignore_folders=[RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR])
        #logging.info("making average from following volumes\n {}".format('\n'.join(vols)))

        average = common.average(vols)

        sitk.WriteImage(average, out_path, True)


class TargetBasedRegistration(ElastixRegistration):
    def __init__(self, *args):
        super(TargetBasedRegistration, self).__init__(*args)
        self.fixed = None
        self.fix_folding = False

    def set_target(self, target):
        self.fixed = target

    def run(self):

        if self.movdir.is_file():
            moving_imgs = [self.movdir]
        else:
            moving_imgs = common.get_file_paths(self.movdir, ignore_folders=[RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR])  # This breaks if not ran from config dir

        if len(moving_imgs) < 1:
            raise common.LamaDataException("No volumes in {}".format(self.movdir))

        for mov in moving_imgs:
            mov_basename = mov.stem
            outdir = self.stagedir / mov_basename
            outdir.mkdir(parents=True)

            cmd = {'mov': str(mov),
                   'fixed': str(self.fixed),
                   'outdir': str(outdir),
                   'elxparam_file': str(self.elxparam_file),
                   'threads': self.threads,
                   'fixed': str(self.fixed)}
            if self.fixed_mask is not None:
                cmd['fixed_mask'] = str(self.fixed_mask)

            run_elastix(cmd)

            # Rename the registered output.
            if self.rename_output:
                elx_outfile = outdir / f'result.0.{self.filetype}'
                new_out_name = outdir / f'{mov_basename}.{self.filetype}'

                try:
                    shutil.move(elx_outfile, new_out_name)
                except IOError:
                    logging.error('Cannot find elastix output. Ensure the following is not set: (WriteResultImage  "false")')
                    raise

                move_intemediate_volumes(outdir)

            # add registration metadata
            reg_metadata_path = outdir / common.INDV_REG_METADATA
            fixed_vol_relative = relpath(self.fixed, outdir)
            reg_metadata = {'fixed_vol': fixed_vol_relative}

            with open(reg_metadata_path, 'w') as fh:
                fh.write(yaml.dump(reg_metadata, default_flow_style=False))

            if self.fix_folding:
                # Remove any folds folds in the Bsplines, overwtite inplace
                tform_param_file = outdir / ELX_TRANSFORM_NAME
                unfold_bsplines(tform_param_file, tform_param_file)

                # Retransform the moving image with corrected tform file
                cmd = [
                    'transformix',
                    '-in', str(mov),
                    '-out', str(outdir),
                    '-tp', tform_param_file
                ]
                subprocess.call(cmd)
                unfolded_moving_img = outdir / 'result.nrrd'
                new_out_name.unlink()
                shutil.move(unfolded_moving_img, new_out_name)


class PairwiseBasedRegistration(ElastixRegistration):

    def __init__(self, *args):
        super(PairwiseBasedRegistration, self).__init__(*args)
        self.inputs_and_mean_tp = {}

    def run(self):

        # If inputs_vols is a file get the specified root and paths from it
        if isdir(self.movdir):
            movlist = common.get_file_paths(self.movdir)
        else:
            movlist = common.get_inputs_from_file_list(self.movdir, self.config_dir)

        if len(movlist) < 1:
            raise common.LamaDataException("No volumes in {}".format(self.movdir))

        for fixed in movlist:  # Todo: change variable name fixed to moving
            tp_file_paths = defaultdict(list)
            full_tp_file_paths = []
            fixed_basename = splitext(basename(fixed))[0]
            fixed_dir = self.paths.make(join(self.stagedir, fixed_basename), 'f')

            for moving in movlist:
                if basename(fixed) == basename(moving):
                    continue
                moving_basename = splitext(basename(moving))[0]
                outdir = join(fixed_dir, moving_basename)
                common.mkdir_force(outdir)

                run_elastix({'mov': moving,
                             'fixed': fixed,
                             'outdir': outdir,
                             'elxparam_file': self.elxparam_file,
                             'threads': self.threads,
                             'fixed': fixed})
                # Get the resolution tforms
                tforms = list(sorted([x for x in os.listdir(outdir) if x .startswith(RESOLUTION_TP_PREFIX)]))
                # get the full tform that spans all resolutions
                full_tp_file_paths.append(join(outdir, FULL_STAGE_TP_FILENAME))

                # Add the tforms to a resolution-specific list so we can generate deformations from any range
                # of deformations later
                for i, tform in enumerate(tforms):
                    tp_file_paths[i].append(join(outdir, tform))

                # add registration metadata
                reg_metadata_path = join(outdir, common.INDV_REG_METADATA)
                fixed_vol_relative = relpath(fixed, outdir)
                reg_metadata = {'fixed_vol': fixed_vol_relative}
                with open(reg_metadata_path, 'w') as fh:
                    fh.write(yaml.dump(reg_metadata, default_flow_style=False))

            for i, files_ in tp_file_paths.items():
                mean_tfom_name = "{}{}.txt".format(RESOLUTION_TP_PREFIX, i)
                self.generate_mean_tranform(files_, fixed, fixed_dir, mean_tfom_name, self.filetype)
            self.generate_mean_tranform(full_tp_file_paths, fixed, fixed_dir, FULL_STAGE_TP_FILENAME, self.filetype)

    @staticmethod
    def generate_mean_tranform(tp_files, fixed_vol, out_dir, tp_out_name, filetype):
        """

        Parameters
        ----------
        tp_files: dict
            each entry contains a list of tforms from a specifc resolution
            eg: {0: [tpfile, tpfile]}
        fixed_vol: str
            path to fixed volume
        out_dir: str
            path to output directory
        filetype
            Filetype for output
        -------

        """
        # get the first tp file to use as template
        template = tp_files[0]
        mean_tp_file = join(out_dir, tp_out_name)

        with open(template, 'r') as tf, open(mean_tp_file, 'w') as outf:
            for line in tf:
                if line.startswith('(Transform '):
                    line = '(Transform "WeightedCombinationTransform")\n'
                if line.startswith('(NumberOfParameters'):
                    line = '(NumberOfParameters {})\n'.format(len(tp_files))
                elif line.startswith('(TransformParameters'):
                    tfparms_str = str('1 ' * len(tp_files)).strip()
                    line = '(TransformParameters {})\n'.format(tfparms_str)
                outf.write(line)
            outf.write('(SubTransforms {})\n'.format(' '.join('"{0}"'.format(x) for x in tp_files)))
            outf.write('(NormalizeCombinationWeights "true")\n')

        cmd = ['transformix',
               '-in', fixed_vol,
               '-tp', mean_tp_file,
               '-out', out_dir,
               ]
        try:
            subprocess.check_output(cmd)
        except Exception as e:  # Can't seem to log CalledProcessError
            logging.warn('transformix failed {}'.format(', '.join(cmd)))
            raise RuntimeError('### Transformix failed creating average ###\nelastix command:{}'.format(cmd))
        else:
            # rename the image output
            bname = splitext(basename(fixed_vol))[0]
            elx_outfile = join(out_dir, 'result.' + filetype)
            new_out_name = join(out_dir, '{}.{}'.format(bname, filetype))
            shutil.move(elx_outfile, new_out_name)


def run_elastix(args):
    cmd = ['elastix',
           '-f', args['fixed'],
           '-m', args["mov"],
           '-out', args['outdir'],
           '-p', args['elxparam_file'],
           ]

    if args.get('threads'):
        cmd.extend(['-threads', str(args['threads'])])

    if args.get('fixed_mask'):
        cmd.extend(['-fMask', args['fixed_mask']])

    try:
        a = subprocess.check_output(cmd)
    except Exception as e:  # can't seem to log CalledProcessError
        logging.exception('registration falied:\n\ncommand: {}\n\n error:{}'.format(cmd, e.output))
        raise


def move_intemediate_volumes(reg_outdir: Path):
    """
    If using elastix multi-resolution registration and outputting image each resolution, put the intermediate files
    in a separate folder.

    Do the same with pyramid images
    """

    intermediate_imgs = list(reg_outdir.rglob('*result.0.R*'))  #[x for x in imgs if basename(x).startswith('result.')]
    if len(intermediate_imgs) > 0:

        reolution_img_dir = reg_outdir / RESOLUTION_IMGS_DIR
        common.mkdir_force(reolution_img_dir)
        for int_img in intermediate_imgs:
            shutil.move(str(int_img), str(reolution_img_dir))

    pyramid_imgs = list(reg_outdir.rglob('*ImagePyramid*'))
    if len(pyramid_imgs) > 0:

        img_pyramid_dir = reg_outdir / IMG_PYRAMID_DIR
        common.mkdir_force(img_pyramid_dir)
        for pyr_img in pyramid_imgs:
            shutil.move(str(pyr_img), str(img_pyramid_dir))



