"""
13/09/19

Display specimens from a line along with stage-matched baselines
"""

import sys
from vpv.utils import data_loader_2
import pandas as pd
from pathlib import Path
from lama.common import STAGING_INFO_FILENAME

wt_dir = Path('/mnt/IMPC_research/neil/E14.5/baselines/output')
mut_dir = Path('/mnt/IMPC_research/neil/E14.5/mutants/output')

baseline_df = pd.read_csv(wt_dir / STAGING_INFO_FILENAME)
mutant_df = pd.read_csv(mut_dir / STAGING_INFO_FILENAME)


line = 'nxn'


lin_specs = mut_dir[mut_dir.line == line]

if len(lin_specs) < 1:
    print(f'No specimens available for {line}')

