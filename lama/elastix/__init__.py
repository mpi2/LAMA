from lama import common
# common.add_elastix_env()
# common.test_installation('transformix')

ELX_TRANSFORM_PREFIX = 'TransformParameters.0.txt'
ELX_PARAM_PREFIX = 'elastix_params_'
ELX_INVERTED_POINTS_NAME = 'outputpoints.vtk'
FILE_FORMAT = '.nrrd'
LOG_FILE = 'inversion.log'
TRANSFORMIX_OUT_NAME = 'result.nrrd'
INVERSION_DIR_NAME = 'Inverted_transform_parameters'
LABEL_INVERTED_TRANFORM = 'labelInvertedTransform.txt'
IMAGE_INVERTED_TRANSFORM = 'ImageInvertedTransform.txt'
VOLUME_CALCULATIONS_FILENAME = "organvolumes.csv"
INVERT_CONFIG = 'invert.yaml'
IGNORE_FOLDER = 'resolution_images'  # When reading images from dir and subdirs, ignore images in this folder
TRANSFORMIX_OUT = 'result.nrrd'  # This will not be correct if filetype is not nrrd
