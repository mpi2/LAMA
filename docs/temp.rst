=====================================================
MRC Harwel Registration pipeline
=====================================================

mrch_regpipeline.py

:Author:
  `Neil Horner`

:Organization:
  Medical Research Council (MRC) Harwell, Oxfordshire, UK

:Version: 1.0.0

This is the main module of the registration of the mouse embryo registration pipeline.  It can be used for creating population averages and for the detection of anatomical phenotypes ( see pheno_detect.py)



Requirements
------------
* Python 2.7
* PyYAML
* SimpleITK

Usage
-----

.. code-block:: bash

    ./mrch_regpipeline.py -c wildtype_config.yaml


