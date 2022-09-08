#!/usr/bin/env bash

script_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

pytest -m "not notest" "$script_dir"/test_data_generation.py; # Run this first as the stats tests need the output
pytest -m "not notest" "$script_dir"/test_standard_stats.py;
pytest -m "not notest" "$script_dir"/test_permutation_stats.py;