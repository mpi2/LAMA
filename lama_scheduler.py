#!/usr/bin/env python

"""

"""

import time
import lama
import pheno_detect


def lama_job_runner(job_file, freq=10):
    while True:
        try:
            with open(job_file) as fh:
                all_jobs = fh.readlines()
                if not all_jobs:  # no jobs in file
                    time.sleep(freq)
                    continue
                current_job = all_jobs[0].strip().split(' ')
                command = current_job[0]
                args = current_job[1:]
        except IOError:
            print 'File may be open'
            time.sleep(freq)
            continue  # The file may be open for reading so try in a bit
        if command == 'lama':
            # remove current job from list
            del all_jobs[0]
            write_remaining_jobs(all_jobs, job_file)
            lama.RegistraionPipeline(args[0])
        if command == 'phenodetect':
            # remove current job from list
            del all_jobs[0]
            write_remaining_jobs(all_jobs, job_file)
            pheno_detect.PhenoDetect(*args)


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
    parser.add_argument('-f', '--freq', dest='freq', help='how often to check for new jobs (seconds)', required=False,
                        type=int, default=10)
    args = parser.parse_args()
    lama_job_runner(args.job_file, args.freq)
