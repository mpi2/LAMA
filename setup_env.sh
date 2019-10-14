#!/usr/bin/env bash

pip3 install --user pipenv
pipenv --python 3.6 #

# Install python dependencies from lama/PipFile
pipenv install
