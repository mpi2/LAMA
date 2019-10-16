# coding: utf-8

from setuptools import setup, find_packages
from pathlib import Path

setup(
    name='lama_phenotype_detection',
    download_url='https://github.com/mpi2/lama/archive/0.9.tar.gz',
    version='0.9.1',
    packages=find_packages(exclude=("dev",)),
    python_requires='>=3.6.*',
    install_requires=[
        'appdirs',
        'matplotlib>=2.2.0',
        'numpy>=1.15.0',
        'pandas>=0.23.4',
        'scikit-learn>=0.19.2',
        'scipy>=1.1.0',
        'scikit-image>=0.15.0',
        'seaborn>=0.9.0',
        'statsmodels>=0.9.0',
        'PyYAML>=3.13',
        'SimpleITK>=1.1.0',
        'filelock',
        'psutil',
        'logzero',
        'addict',
        'toml',
        'pynrrd',
        'pytest'
    ],
    extras_require={
        'dev': ['pyradiomics'],
    },
    url='https://www.har.mrc.ac.uk/',
    license='Apache2',
    author='Neil Horner',
    author_email='n.horner@har.mrc.ac.uk, bit@har.mrc.ac.uk',
    description='Phenotype detection pipeline for finding abnormalities in mouse embryos',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
     ],
    keywords=['image processing', 'bioinformatics', 'phenotype'],
    entry_points ={
            'console_scripts': [
                'lama_reg=scripts.lama_reg:main',
                'lama_get_test_data=scripts.lama_get_test_data:main',
                'lama_get_tutorial_data=scripts.lama_get_tutorial_data:main',
                'lama_job_runner=scripts.lama_job_runner:main',
                'lama_permutation_stats=scripts.lama_permutation_stats:main',
                'lama_stats=scripts.lama_stats:main'
            ]
        },
)
