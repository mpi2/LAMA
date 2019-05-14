"""
Run some tests on the lama config validator
"""

from pathlib import Path
from lama.registration_pipeline.validate_config import LamaConfig, LamaConfigError

test_data_root = join(current_dir, 'test_data')
registration_root = join(test_data_root, 'registration_test_data')
test_config_path  = Path(registration_root) / 'registration_config.toml'