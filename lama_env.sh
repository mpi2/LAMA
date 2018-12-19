#! /bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
SCRIPTS="/scripts"
UTILS="/utilities"
script_dir=$DIR$SCRIPTS
utils_dir=$DIR$UTILS

# Add the script directories to the PATH
export PATH=$utils_dir:$script_dir:$PATH
# Actiavte the python virtual environment
pipenv shell



