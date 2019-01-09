#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"


FN="/../lama/stats/standard_stats/lama_stats_new.py"

file=$LAMA_SCRIPT_DIR$FN
echo $file
python3 $file