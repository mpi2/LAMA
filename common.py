import logging
import subprocess
import os
import datetime
import shutil

LOG_MODE = logging.DEBUG

def init_log(logpath, name):

    logging.basicConfig(filename=logpath, level=LOG_MODE, filemode="w")
    logging.info(name)

    # switch working dir to this module's location
    orig_wd = os.getcwd()
    module_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(module_dir)

    try:
        git_log = subprocess.check_output(['git', 'log', '-n', '1'])
        git_commit = git_log.splitlines()[0]
        git_branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
    except subprocess.CalledProcessError:
        git_commit = "Git commit info not available"
        git_branch = "Git branch not available"

    logging.info("git branch: {}".format(git_branch.strip()))
    logging.info("git {}".format(git_commit.strip()))

    # back to original working dir
    os.chdir(orig_wd)

def log_time(msg):
    now = datetime.datetime.now()
    logging.info("{}: {}/{}/{} - {}:{}".format(msg, now.day, now.month, now.year, now.hour, now.minute))

def mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)