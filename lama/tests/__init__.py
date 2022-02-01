import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

# When running lama,. first set the environment up using 'source lama_env.sh'
# But for the tests, we shall just set it here
import sys
import yaml

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
stats_config_dir = current_dir / 'configs'
permutation_stats_dir = stats_test_data_dir / 'permutation_stats'

qc_flags_dir = registration_root / 'qc_flag_files'


# current_dir = Path.cwd()


INPUT_DIR = test_data_root / 'input_data'

baseline_input_dir = registration_root / 'baseline'
population_test_dir = test_data_root / 'population_average_data'


#for Kyle's perm testing

cfg_path = Path(
    "C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/LAMA/lama/tests/configs/permutation_stats/perm_no_qc.yaml")
def p(path):
    if path is None:
        return

    cfg_dir = Path(cfg_path).parent

    resolved = (cfg_dir / path).resolve()

    print(resolved)
    if not resolved.exists():
        raise FileNotFoundError(f'Cannot find: {resolved}')
    return resolved


with open(cfg_path, 'r') as fh:
    cfg = yaml.load(fh)

wt_dir = p(cfg['wildtype_dir'])
mut_dir = p(cfg['mutant_dir'])

out_dir = p(cfg.get('output_dir', Path(cfg_path).parent))



n_perm = int(cfg.get('n_permutations', 1000))
label_meta = p(cfg.get('label_metadata'))
label_map = p(cfg.get('label_map'))
wev_norm = bool(cfg.get('norm_to_whole_embryo_vol', True))
qc_file = p(cfg.get('qc_file'))
voxel_size = float(cfg.get('voxel_size', 1.0))

treat_dir = p(cfg['treatment_dir'])
inter_dir = p(cfg['interaction_dir'])
two_way = bool(cfg.get('two_way', False))

