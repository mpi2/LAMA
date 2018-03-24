from os.path import join, dirname, realpath

current_dir = dirname(realpath(__file__))
test_data_root = join(current_dir, '..', 'test_data', 'registration_test_data')

INPUT_DIR = join(test_data_root, 'input_data')