#!/bin/bash

../../mrch_regpipeline.py -c wt/wt_config.yaml && ../../pheno_detect.py -c wt/wt_config.yaml -i mut/inputs -p mut
