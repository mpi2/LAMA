#!/usr/bin/env bash
 script_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
 pytest -m "not notest" $script_dir