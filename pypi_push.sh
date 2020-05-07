#! /bin/sh

python3.6 setup.py sdist bdist_wheel 
# Need to set repository info in ~/.pypirc
twine upload dist/*
