from os.path import join
import os
import common
import logging
import sys

def validate_config(self, config):
        """
        Do some checks on the config file to check for errors
        :param config:
        :return: list, empty or containing errors


        TODO: Add stats checking. eg. if using lmR, need to specify formula
        """
        report = []
        required_params = ['output_dir', 'fixed_volume',
                           'inputvolumes_dir', 'filetype',
                           'global_elastix_params',
                           'registration_stage_params']

        for p in required_params:
            if p not in config:
                report.append("Entry '{}' is required in the config file".format(p))

        stages = config['registration_stage_params']

        if len(stages) < 1:
            report.append("No stages specified")

        # Check whether images are 16 bit and if so whether internal representation is set to float
        img_dir = join(self.config_dir, config['inputvolumes_dir'])
        imgs = os.listdir(img_dir)

        logging.info('validating input volumes')
        for im_name in imgs:
            image_path = join(img_dir, im_name)
            if not os.path.isfile(image_path):
                logging.info('Something wrong with the inputs. Cannot find {}'.format(image_path))
            array = common.img_path_to_array(image_path)
            if array.dtype in (np.int16, np.uint16):
                try:
                    internal_fixed = config['global_elastix_params']['FixedInternalImagePixelType']
                    internal_mov = config['global_elastix_params']['MovingInternalImagePixelType']
                    if internal_fixed != 'float' or internal_mov != 'float':
                        raise TypeError
                except (TypeError, KeyError):
                    logging.error("If using 16 bit input volumes, 'FixedInternalImagePixelType' and 'MovingInternalImagePixelType should'" \
                                  "be set to 'float' in the global_elastix_params secion of the config file")
                    sys.exit()

        if len(report) > 0:
            for r in report:
                logging.error(r)
            sys.exit()

