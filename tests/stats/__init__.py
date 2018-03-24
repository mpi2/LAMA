from os.path import join, dirname, realpath

current_dir = dirname(realpath(__file__))
test_data_root = join(current_dir, '..', 'test_data', 'stats_test_data')

CONFIG_DIR = join(test_data_root, 'config_files')
INPUT_DIR = join(test_data_root, 'input_data')
