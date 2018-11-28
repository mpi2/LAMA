import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

current_dir = Path(__file__).parent
test_data_root = current_dir / 'test_data'
registration_data_root = current_dir / 'test_data' /'registration_test_data'

stats_test_data_dir = test_data_root / 'stats_test_data'
stats_config_dir = stats_test_data_dir / 'config_files'