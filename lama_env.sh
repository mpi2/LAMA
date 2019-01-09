#! /bin/bash

# If running on Linux, setup environment to add aliases variables so lamam scripts can be run from anywhere

# Get the root directory of the project
LAMA_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# Add the root directory to the PYTHONPATH so the LAMA pkagae can be imported regardless of where calling from.
export PYTHONPATH="$LAMA_ROOT:$PYTHONPATH"

SCRIPTS="/scripts/"
UTILS="/utilities/"

LAMA_SCRIPT_DIR=$LAMA_ROOT$SCRIPTS
LAMA_UTILS_DIR=$LAMA_ROOT$UTILS

echo "Setting lama script environment variables"

INT="python3 "
PYTHON_FILE=$LAMA_SCRIPT_DIR$NAME

EX_PREFIX="pipenv run $INT$LAMA_SCRIPT_DIR"

alias lama_run="${EX_PREFIX}lama_reg.py"
alias lama_stats="${EX_PREFIX}lama_stats.py"
alias lama_job_runner="${EX_PREFIX}job_runner.py"
alias lama_get_test_data="${EX_PREFIX}lama_get_test_data.py"
alias lama_get_tutorial_data="${EX_PREFIX}lama_get_tutorial_data.py"



