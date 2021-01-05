import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

# When running lama,. first set the environment up using 'source lama_env.sh'
# But for the tests, we shall just set it here
import sys
lama_dir = Path.cwd().parent
sys.path.insert(0, lama_dir)

current_dir = Path(__file__).parent
test_data_root = current_dir / 'test_data'
registration_root = current_dir / 'test_data' / 'mutant_and_baseline_data'
wt_registration_dir = registration_root / 'baseline'
mut_registration_dir = registration_root / 'mutant'
target_dir = registration_root / 'target'
target_img = target_dir / '301015_deformable_to_8_rescaled_8bit.nrrd'

stats_test_data_dir = test_data_root / 'stats_test_data'
stats_output_dir = stats_test_data_dir / 'test_output'
stats_config_dir = stats_test_data_dir / 'config_files'
permutation_stats_dir = stats_test_data_dir / 'permutation_stats'

qc_flags_dir = registration_root / 'qc_flag_files'


# current_dir = Path.cwd()


INPUT_DIR = test_data_root / 'input_data'

baseline_input_dir = registration_root / 'baseline'
population_test_dir = test_data_root / 'population_average_data'



