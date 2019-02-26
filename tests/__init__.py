import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

# When running lama,. first set the environment up using 'spurce lama_env.sh'
# But for the tests, we shall just set it here
import sys
lama_dir = Path.cwd().parent
sys.path.insert(0, lama_dir)

current_dir = Path(__file__).parent
test_data_root = current_dir / 'test_data'
registration_data_root = current_dir / 'test_data' /'registration_test_data'
wt_registration_dir = registration_data_root / 'baseline'
mut_registration_dir = registration_data_root / 'mutant'
target_dir = registration_data_root / 'target'

stats_test_data_dir = test_data_root / 'stats_test_data'
stats_output_dir = stats_test_data_dir / 'test_output'
stats_config_dir = stats_test_data_dir / 'config_files'


