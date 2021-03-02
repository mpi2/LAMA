#!/usr/bin/env bash
 script_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
 pytest -q -m "not notest" $script_path