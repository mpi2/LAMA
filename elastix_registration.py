import common
import logging
from os.path import join, isdir, splitext, basename, relpath, exists
import subprocess
from multiprocessing import Pool
import shutil
import yaml
import SimpleITK as sitk

INDV_REG_METADATA = 'reg_metadata.yaml'
TP_FILENAME = 'TransformParameters.0.txt'


#TODO: Need to upadte pairwise to accunt for changes in threading

class ElastixRegistration(object):

    def __init__(self, elxparam_file, movdir, stagedir, paths, filetype, threads=None, fixed_mask=None,
                 lama_multithread=False):
        self.elxparam_file = elxparam_file
        self.movdir = movdir
        self.stagedir = stagedir
        self.fixed_mask = fixed_mask
        self.paths = paths
        self.filetype = filetype
        self.threads = int(threads)
        self.lama_multithread = lama_multithread

    def make_average(self, out_path):
        """
        Create an average of the the input embryo volumes.
        This will search subfolders for all the registered volumes within them
        """
        vols = common.GetFilePaths(self.stagedir)

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
        if self.lama_multithread:
            self.run_multithread()
        else:
            self.run_single_thread()

    def run_single_thread(self):

        if self.fixed is None:
            raise NameError('In TargetBasedRegistration, target must be set using set_target ')
        # If inputs_vols is a file get the specified root and paths from it
        if isdir(self.movdir):
            movlist = common.GetFilePaths(self.movdir)
        else:
            movlist = common.get_inputs_from_file_list(self.movdir)

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
            shutil.move(elx_outfile, new_out_name)

            # add registration metadata
            reg_metadata_path = join(outdir, INDV_REG_METADATA)
            fixed_vol_relative = relpath(self.fixed, outdir)
            reg_metadata = {'fixed_vol': fixed_vol_relative}
            with open(reg_metadata_path, 'w') as fh:
                fh.write(yaml.dump(reg_metadata, default_flow_style=False))

    def run_multithread(self):
        if self.fixed is None:
            raise NameError('In TargetBasedRegistration, target must be set using set_target ')
        # If inputs_vols is a file get the specified root and paths from it
        if isdir(self.movdir):
            movlist = common.GetFilePaths(self.movdir)
        else:
            movlist = common.get_inputs_from_file_list(self.movdir)

        jobs = []
        for mov in movlist:
            mov_basename = splitext(basename(mov))[0]
            outdir = self.paths.make(join(self.stagedir, mov_basename), 'f')

            job = {'mov': mov,
                   'fixed': self.fixed_mask,
                   'outdir': outdir,
                   'elxparam_file': self.elxparam_file,
                   'threads': self.threads,
                   'fixed': self.fixed,
                   'fixed_mask': self.fixed_mask,
                   'mov_base': mov_basename}
            jobs.append(job)

        pool = Pool(self.threads)
        try:
            pool.map(run_elastix, jobs)
        except KeyboardInterrupt:
            print 'terminating inversion'
            pool.terminate()
            pool.join()

        # now rename the registered outputs
        self.rename_multi_output(jobs)

    def rename_multi_output(self, jobs):
        # Rename the registered output.
        for job in jobs:
            outdir = job['outdir']
            mov_basename = job['mov_base']
            elx_outfile = join(outdir, 'result.0.{}'.format(self.filetype))
            new_out_name = join(outdir, '{}.{}'.format(mov_basename, self.filetype))
            shutil.move(elx_outfile, new_out_name)

            # add registration metadata
            reg_metadata_path = join(outdir, INDV_REG_METADATA)
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
            movlist = common.get_inputs_from_file_list(self.movdir)

        for fixed in movlist:  # Todo: change variable name fixed to moving
            tp_file_paths = []
            fixed_basename = splitext(basename(fixed))[0]
            fixed_dir = self.paths.make(join(self.stagedir, fixed_basename), 'f')

            for moving in movlist:
                if basename(fixed) == basename(moving):
                    continue
                moving_basename = splitext(basename(moving))[0]
                outdir = join(fixed_dir, moving_basename)
                common.mkdir_force(outdir)

                self.run_elastix(fixed, moving, outdir)  #  Flipped the moving and fixed to see if we can get around inverting transforms
                tp_file_paths.append(join(outdir, TP_FILENAME))

                # add registration metadata
                reg_metadata_path = join(outdir, INDV_REG_METADATA)
                fixed_vol_relative = relpath(fixed, outdir)
                reg_metadata = {'fixed_vol': fixed_vol_relative}
                with open(reg_metadata_path, 'w') as fh:
                    fh.write(yaml.dump(reg_metadata, default_flow_style=False))

            mean_tp_file = self.generate_mean_tranform(tp_file_paths, fixed, fixed_dir)

            self.inputs_and_mean_tp[fixed] = mean_tp_file

    @staticmethod
    def generate_mean_tranform(tp_files, input_vol, out_dir):
        # get the first tp file to use as template
        template = tp_files[0]
        mean_tp_file = join(out_dir, 'meanTransformParameter.txt')
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
               '-in', input_vol,
               '-tp', mean_tp_file,
               '-out', out_dir,
               ]
        try:
            subprocess.check_output(cmd)
        except Exception as e:  # Can't seem to log CalledProcessError
            logging.warn('transformix failed {}'.format(', '.join(cmd)))
            raise RuntimeError('### Transformix failed creating average ###\nelastix command:{}'.format(cmd))
        else:
            bname = splitext(basename(input_vol))[0]
            elx_outfile = join(out_dir, 'result.nrrd')
            new_out_name = join(out_dir, '{}.nrrd'.format(bname))
            shutil.move(elx_outfile, new_out_name)


def run_elastix(args):
    cmd = ['elastix',
           '-f', args['fixed'],
           '-m', args['mov'],
           '-out', args['outdir'],
           '-p', args['elxparam_file'],
           ]

    if args['threads']:
        cmd.extend(['-threads', str(args['threads'])])

    if args['fixed_mask']:
        cmd.extend(['-fMask', args['fixed_mask']])

    try:
        subprocess.check_output(cmd)
    except Exception as e:  # can't seem to log CalledProcessError
        logging.exception('registration falied:\n\ncommand: {}\n\n error:{}'.format(cmd, e))
        raise