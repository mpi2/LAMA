import common
import logging
import sys
from os.path import join, isdir, splitext, basename, relpath, exists
import subprocess
import os
import shutil
import yaml
import SimpleITK as sitk
from utilities import batch_convert_images
from collections import defaultdict


REOLSUTION_TP_PREFIX = 'TransformParameters.0.R'
FULL_STAGE_TP_FILENAME = 'TransformParameters.0.txt'
RESOLUTION_IMG_FOLDER = 'resolution_images'



#TODO: Need to upadte pairwise to accunt for changes in threading

class ElastixRegistration(object):

    def __init__(self, elxparam_file, movdir, stagedir, paths, filetype, threads=None, fixed_mask=None, config_dir=None):
        self.elxparam_file = elxparam_file
        self.movdir = movdir
        self.stagedir = stagedir
        self.fixed_mask = fixed_mask
        self.paths = paths
        self.filetype = filetype
        self.threads = int(threads)
        self.config_dir = config_dir
        # A subset of volumes from folder to register

    def make_average(self, out_path):
        """
        Create an average of the the input embryo volumes.
        This will search subfolders for all the registered volumes within them
        """
        vols = common.GetFilePaths(self.stagedir, ignore_folder=RESOLUTION_IMG_FOLDER)
        logging.info("making average from following volumes\n {}".format('\n'.join(vols)))

        average = common.Average(vols)

        sitk.WriteImage(average, out_path, True)  # Compressed=True
        # Check that it's been created
        if not exists(out_path):
            logging.error('Cannot make average at {}'.format(out_path))


class TargetBasedRegistration(ElastixRegistration):
    def __init__(self, *args):
        super(TargetBasedRegistration, self).__init__(*args)
        self.fixed = None

    def set_target(self, target):
        self.fixed = target

    def run(self):
        self.run_single_thread()

    def run_single_thread(self):

        # If inputs_vols is a file get the specified root and paths from it
        if isdir(self.movdir):
            movlist = common.GetFilePaths(self.movdir, ignore_folder=RESOLUTION_IMG_FOLDER)  # This breaks if not ran from config dir
        else:
            movlist = common.get_inputs_from_file_list(self.movdir, self.config_dir)

        for mov in movlist:
            mov_basename = splitext(basename(mov))[0]
            outdir = self.paths.make(join(self.stagedir, mov_basename), 'f')

            run_elastix({'mov': mov,
                         'fixed': self.fixed,
                         'outdir': outdir,
                         'elxparam_file': self.elxparam_file,
                         'threads': self.threads,
                         'fixed': self.fixed,
                         'fixed_mask': self.fixed_mask})

            # Rename the registered output.
            elx_outfile = join(outdir, 'result.0.{}'.format(self.filetype))
            new_out_name = join(outdir, '{}.{}'.format(mov_basename, self.filetype))
            try:
                shutil.move(elx_outfile, new_out_name)
            except IOError:
                logging.error('Cannot find elastix output. Is the following set (WriteResultImage  "false")')
                sys.exit(1)
            move_intemediate_volumes(outdir)

            # add registration metadata
            reg_metadata_path = join(outdir, common.INDV_REG_METADATA)
            fixed_vol_relative = relpath(self.fixed, outdir)
            reg_metadata = {'fixed_vol': fixed_vol_relative}
            with open(reg_metadata_path, 'w') as fh:
                fh.write(yaml.dump(reg_metadata, default_flow_style=False))


class PairwiseBasedRegistration(ElastixRegistration):

    def __init__(self, *args):
        super(PairwiseBasedRegistration, self).__init__(*args)
        self.inputs_and_mean_tp = {}

    def run(self):

        # If inputs_vols is a file get the specified root and paths from it
        if isdir(self.movdir):
            movlist = common.GetFilePaths(self.movdir)
        else:
            movlist = common.get_inputs_from_file_list(self.movdir, self.config_dir)

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
                tforms = list(sorted([x for x in os.listdir(outdir) if x .startswith(REOLSUTION_TP_PREFIX)]))
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

            for i, files_ in tp_file_paths.iteritems():
                mean_tfom_name = "{}{}.txt".format(REOLSUTION_TP_PREFIX, i)
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
            #outf.write('(NormalizeCombinationWeights "true")\n')

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
           '-m', args['mov'],
           '-out', args['outdir'],
           '-p', args['elxparam_file'],
           ]

    if args.get('threads'):
        cmd.extend(['-threads', str(args['threads'])])

    if args.get('fixed_mask'):
        cmd.extend(['-fMask', args['fixed_mask']])

    try:
        subprocess.check_output(cmd)
    except Exception as e:  # can't seem to log CalledProcessError
        logging.exception('registration falied:\n\ncommand: {}\n\n error:{}'.format(cmd, e))
        raise


def move_intemediate_volumes(reg_outdir):
    """
    If using elastix multi-resolution registration and outputing image each resolution, put the intemediate files
    in a
    Returns
    -------
    """
    imgs = common.GetFilePaths(reg_outdir)
    intermediate_imgs = [x for x in imgs if basename(x).startswith('result.')]
    if len(intermediate_imgs) > 0:
        int_dir = join(reg_outdir, RESOLUTION_IMG_FOLDER)
        common.mkdir_force(int_dir)
        for int_img in intermediate_imgs:
            shutil.move(int_img, int_dir)
        batch_convert_images.cast_and_rescale_to_8bit(int_dir, int_dir)

