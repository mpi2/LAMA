"""
lama_jobrunner.py allows us to run a bunch of jobs, but this has to be run on multiple machine.
If all we want to do is farm the jobs out across a griod, we can use this script. Only run this on one machine and it
will call jobrunner multuple times.
"""

import toml
import sys
import fabric


def main():
    """
    This function exists so that arguments can be taken from sys.argv if running as a setuptools entry point script.
    """
    config_path = sys.argv[1]
    run_on_grid(config_path)


def run_on_grid(config_path):

    cfg = toml.load(config_path)
    cmd = f'{cfg["grid_cmd"]} "{cfg["docker_cmd"]} \'{cfg["lama_cmd"]}\'"'


if __name__ == '__main__':
    config_path = sys.argv[1]
    run_on_grid(config_path)
    run_on_grid(config_path)

