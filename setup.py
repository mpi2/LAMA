# coding: utf-8

from setuptools import setup, find_packages

setup(
    name='lama_phenotype_detection',
    download_url='https://github.com/mpi2/lama/archive/0.9.70.tar.gz',
    version='0.9.70',
    packages=find_packages(exclude=("dev")),
    package_data={'': ['current_commit',
                       'stats/rscripts/lmFast.R',
                       'stats/rscripts/r_padjust.R']},  # Puts it in the wheel dist. MANIFEST.in gets it in source dist
    include_package_data=True,
    install_requires=[
        'appdirs',
        'matplotlib>=2.2.0',
        'numpy>=1.15.0',
        'pandas>=1.1.0',
        'scikit-learn>=0.19.2',
        'scipy>=1.1.0',
        'scikit-image==0.17.2',
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
    url='https://github.com/mpi2/LAMA',
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
                'lama_reg=lama.scripts.lama_reg:main',
                'lama_get_test_data=lama.scripts.lama_get_test_data:main',
                'lama_get_walkthrough_data=lama.scripts.lama_get_walkthrough_data:main',
                'lama_job_runner=lama.scripts.lama_job_runner:main',
                'lama_permutation_stats=lama.scripts.lama_permutation_stats:main',
                'lama_stats=lama.scripts.lama_stats:main',
                'lama_pad_volumes=lama.utilities.lama_pad_volumes:main',
                'lama_convert_16_to_8=lama.utilities.lama_convert_16_to_8:main',
                'lama_img_info=lama.utilities.lama_img_info:main'
            ]
        },
)
