"""
This entry point to run lama stats pipeline based on running a single lne againt a series of baselines.
Call elements of the pipeline from here and keep it simple so as to act as an overall vew of the step ij the pipeline
"""

from . import parse_and_validate_config, get_data
from . parse_and_validate_config import LamaStatsData


def run_all(config_path):
    stats_objs = parse_and_validate_config.get_stats_objects()

    for s in stats_objs:
        run_stats_pipeline(s)

def run_stats_pipeline(stats_data: LamaStatsData):

    # setup logging

    # Read and validate config. Return an object that will have all data needed to run a stats analysis
    # This object will also eventaully contain a stats return object


    # get_data
    get_data(stats_data)

    # Add stats method
    add_stats_method



if __name__ == '__main__':
    run_stats_pipeline()