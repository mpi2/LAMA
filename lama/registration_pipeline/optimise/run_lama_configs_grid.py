"""
Given a input folder of lines (could be centres etc) and a folder of lama configs:
* Run each specimen on a grid node with each config file (same config as used for grid job runner
*

Example toml config
-------------------
grid_cmd = 'source /NGS/grid/dist/GE2011.11p1/informatics/common/settings.sh; /grid/dist/GE2011.11p1/bin/linux-x64/qsub -N neil_lama -wd /grid/output -j y -b yes -P SIG -p 0 -o /grid/output -pe nu {}'
docker_cmd = 'docker run --cap-add SYS_ADMIN --cap-add DAC_READ_SEARCH -h=`hostname` cutter:5000/neil_lama:{} bash -ci'
lama_cmd = 'lama_reg -c {}'
root_dir = '/mnt/IMPC_media/LAMA_staging/e15_5/080620_pop_avg/1/test_all_wts_110620'
HOST = 'hampshire'
USER = 'n.horner'
n_thread = '16'
docker_tag = '354'
SGE_root = '/grid/dist/GE2011.11p1'
-------------------

This setup is specific to Harwell infrastructure, but could be easily modified

"""

import fabric
import time
import shutil
from pathlib import Path
import toml


def run(config_path):

    grid_config = toml.load(config_path)
    root = Path(grid_config['root_dir'])
    inputs_dir = root / 'inputs'
    if not inputs_dir.is_dir():
        raise NotADirectoryError

    out_root = root / 'output'
    out_root.mkdir(exist_ok=True)

    config_dir = root / 'configs'
    if not config_dir.is_dir():
        raise NotADirectoryError

    config_done_dir = root / 'configs_done'
    config_done_dir.mkdir(exist_ok=True)

    while True:

        try:
            lama_config_path = list(config_dir.iterdir())[0]
        except IndexError:
            print(f'Waiting for configs')
            time.sleep(2)
            continue

        lama_config = toml.load(lama_config_path)
        shutil.move(lama_config_path, config_done_dir / lama_config_path.name)

        root_reg_dir = out_root / lama_config_path.stem
        root_reg_dir.mkdir(exist_ok=True)

        # Check if multiple dirs exist in inputs (whch maight be different lines/centres for  example
        # If only files, we just run them all together
        if any([x.is_dir() for x in inputs_dir.iterdir()]):
            input_dirs = [x for x in inputs_dir.iterdir() if x.is_dir()]
        else:
            input_dirs = [inputs_dir]

        for line_dir in input_dirs:
            for input_ in line_dir.iterdir():

                print(f'Specimen: {input_.name}, config: {lama_config_path.name}')
                reg_out_dir = root_reg_dir / line_dir.name / input_.stem
                reg_out_dir.mkdir(exist_ok=True, parents=True)

                new_input_dir = reg_out_dir / 'inputs'
                new_input_dir.mkdir(exist_ok=True)
                shutil.copyfile(input_, new_input_dir / input_.name)

                new_config_name = reg_out_dir / lama_config_path.name

                with open(new_config_name, 'w') as fh:
                    toml.dump(lama_config, fh)

                run_on_grid(new_config_name, grid_config)


def run_on_grid(lama_config_path, grid_config):
    # lama_config_path = str(lama_config_path).replace('/mnt', '')
    c = grid_config
    cmd = f'{c["grid_cmd"]} "{c["docker_cmd"]} \'{c["lama_cmd"]}\'"'
    # Now interpolate our values
    cmd = cmd.format(c['n_thread'], c['docker_tag'], lama_config_path)

    conn = fabric.Connection(c['HOST'], user=c['USER'], inline_ssh_env=True)
    conn.run(cmd, env={'SGE_ROOT': '/grid/dist/GE2011.11p1'})
    conn.close()


if __name__ == '__main__':
    import sys
    cfg_path = sys.argv[1]
    run(cfg_path)

