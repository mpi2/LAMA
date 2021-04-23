"""

Usage:  pytest -qq -m "not notest" test_data_generation.py
The use of -m "not notest" is to be able to omit certain tests with the @pytest.mark.notest decorator
"""

from pathlib import Path
from lama.registration_pipeline import run_lama
import logzero

import os
import shutil
import pytest
import logging
from lama.scripts import lama_job_runner
from lama.tests import (registration_root, mut_registration_dir, wt_registration_dir)




