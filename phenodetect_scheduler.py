#!/usr/bin/env python

"""
Given a list of jobs to do, periodically checks this list, runs the jobs and removes it from the list


-j
The jobs file list is a text file with the name of the mutant folder on each line
Each folder should have n 'inputs' subfolder
    adam23
    atp1a2
    daam1
    epas1
    ethe1
    f10
    fzd3
    gfi1b

-c
The config file is the one gerneated by lama when running the baselines and should have '_pheno_detect.yam' suffix
"""

import time
import lama
import pheno_detect
from os.path import abspath


def lama_job_runner(job_file, config_path, freq=10):

    job_file = abspath(job_file)
    config_path = abspath(config_path)

    while True:

        try:
            with open(job_file) as fh:
                all_jobs = fh.readlines()
                if not all_jobs:  # no jobs in file
                    time.sleep(freq)
                    continue
                mutant_folder = all_jobs[0].strip()

        except IOError:
            print 'File may be open'
            time.sleep(freq)  # The file may be open for reading so try in a bit
            continue

        else:
            # remove current job from list
            write_remaining_jobs(all_jobs[1:], job_file)
            try:
                pheno_detect.PhenoDetect(config_path, mutant_folder)
            except BaseException as e:  # sys.exit does not inherit from Exception
                print "Failed phenodetect job: {}\n{}".format(all_jobs[0], str(e))



def write_remaining_jobs(all_jobs, job_file):
        # Write the remaining jobs
        while True:
            time.sleep(10)
            try:
                with open(job_file, 'w') as fh:
                    fh.writelines(all_jobs)
            except IOError:
                print 'Job file may be locked for writing'
                continue  # The file may be open for reading so try in a bit
            else:
                break


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")
    parser.add_argument('-j', '--dir', dest='job_file', help='file_with jobs list watch for new jobs', required=True)
    parser.add_argument('-c', '--config', dest='config', help='_pheno_detect.yaml config file', required=True)
    parser.add_argument('-f', '--freq', dest='freq', help='how often to check for new jobs (seconds)', required=False,
                        type=int, default=10)
    args = parser.parse_args()
    lama_job_runner(args.job_file, args.config, args.freq)
