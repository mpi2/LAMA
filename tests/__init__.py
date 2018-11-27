import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

current_dir = Path(__file__).parent
test_data_root = current_dir / 'test_data'
registration_data_root = current_dir / 'test_data' /'registration_test_data'