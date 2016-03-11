import common
import logging
from os.path import join, isdir, splitext, basename, relpath, exists
import subprocess
import shutil
import yaml
import SimpleITK as sitk

INDV_REG_METADATA = 'reg_metadata.yaml'
TP_FILENAME = 'TransformParameters.0.txt'


class ElastixRegistration(object):

    def __init__(self, elxparam_file, movdir, stagedir, paths, filetype, threads, fixed_mask=None):
        self.elxparam_file = elxparam_file
        self.movdir = movdir
        self.stagedir = stagedir
        self.fixed_mask = fixed_mask
        self.paths = paths
        self.filetype = filetype
        self.threads = threads

    def run_elastix(self, mov, fixed, outdir):
        cmd = ['elastix',
               '-f', fixed,
               '-m', mov,
               '-out', outdir,
               '-p', self.elxparam_file,
               '-threads', self.threads]

        try:
            subprocess.check_output(cmd)
        except Exception as e:  # can't seem to log CalledProcessError
            logging.exception('registration falied:\n\ncommand: {}\n\n error:{}'.format(cmd, e))
            raise

    def make_average(self):
        raise NotImplementedError


class TargetBasedRegistration(ElastixRegistration):
    def __init__(self, *args):
        super(TargetBasedRegistration, self).__init__(*args)
        self.fixed = None

    def set_target(self, target):
        self.fixed = target

    def run(self):

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

            self.run_elastix(mov, self.fixed, outdir)

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

        for fixed in movlist:
            tp_file_paths = []
            fixed_basename = splitext(basename(fixed))[0]
            fixed_dir = self.paths.make(join(self.stagedir, fixed_basename), 'f')

            for moving in movlist:
                if fixed == moving:
                    continue
                moving_basename = splitext(basename(moving))[0]
                outdir = join(fixed_dir, moving_basename)
                common.mkdir_force(outdir)

                self.run_elastix(moving, fixed, outdir)
                tp_file_paths.append(join(outdir, TP_FILENAME))

                # Rename the registered output.
                elx_outfile = join(outdir, 'result.0.{}'.format(self.filetype))
                new_out_name = join(outdir, '{}.{}'.format(moving_basename, self.filetype))
                shutil.move(elx_outfile, new_out_name)

                # add registration metadata
                reg_metadata_path = join(outdir, INDV_REG_METADATA)
                fixed_vol_relative = relpath(fixed, outdir)
                reg_metadata = {'fixed_vol': fixed_vol_relative}
                with open(reg_metadata_path, 'w') as fh:
                    fh.write(yaml.dump(reg_metadata, default_flow_style=False))

            mean_tp_file = self.generate_mean_tranform_tp_file(tp_file_paths, fixed_dir)
            self.inputs_and_mean_tp[fixed] = mean_tp_file

    @staticmethod
    def generate_mean_tranform_tp_file(tp_files, outdir):
        # get the first tp file to use as template
        template = tp_files[0]
        mean_tp_file = join(outdir, 'meanTransformParameter.txt')
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
        return mean_tp_file

    def make_average(self, out_path):

        temp_files = join(splitext(out_path)[0])
        common.mkdir_force(temp_files)
        transformed_vols = []
        for input_vol, mean_tp_file in self.inputs_and_mean_tp.iteritems():
            out_dir = join(temp_files, splitext(basename(input_vol))[0])
            common.mkdir_force(out_dir)
            cmd = ['transformix',
                   '-in', input_vol,
                   '-tp', mean_tp_file,
                   '-out', out_dir,
                   ]
            try:
                subprocess.check_output(cmd)
            except Exception as e:  # Can't seem to log CalledProcessError
                logging.warn('transformix failed {}'.format(', '.join(cmd)))
                raise RuntimeError('### Transformix failed while transforming mask ###\nelastix command:{}'.format(cmd))

            t_vol = join(out_dir, 'result.nrrd')
            transformed_vols.append(t_vol)

        average = common.Average(transformed_vols)

        sitk.WriteImage(average, out_path, True)  # Compressed=True
        # Check that it's been created
        if not exists(out_path):
            logging.error('Cannot make average at {}'.format(out_path))