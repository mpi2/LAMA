# coding: utf-8

from setuptools import setup, find_packages
from pathlib import Path

# Get __verison_dunder without importing lama
version_file = Path(__file__).resolve().parent / 'lama' / 'version.py'
exec(open(version_file).read())


setup(
    name='dorkylever_lama_phenotype_detection',
    download_url=f'https://github.com/dorkylever/LAMA/archive/refs/tags/1.0.1.tar.gz',
    version="1.0.1",
    packages=find_packages(exclude=("dev")),
    package_data={'': ['current_commit',
                       'stats/rscripts/lmFast.R',
                       'stats/rscripts/r_padjust.R']},  # Puts it in the wheel dist. MANIFEST.in gets it in source dist
    include_package_data=True,
    install_requires=[
        'appdirs',
        'setuptools==59.8.0',
        'matplotlib>=2.2.0',
        'numpy==1.21.5',
        'pandas>=1.1.0',
        'scikit-learn==1.0.2',
        'scipy>=1.1.0',
        'scikit-image==0.17.2',
        'seaborn>=0.9.0',
        'statsmodels>=0.9.0',
        'PyYAML>=3.13',
        'catboost==1.1.0',
        'SimpleITK>=2.1.0',
        'pyradiomics>=3.0.1',
        'threadpoolctl==3.1.0',
        'imbalanced-learn==0.9.0',
        'raster-geometry',
        'filelock',
        'psutil==5.9.3',
        'plotly',
        'logzero==1.7.0',
        'addict',
        'toml',
        'pynrrd',
        'pytest',
        'tqdm',
        'gitpython',
        'pacmap',
        'shap',
        'joblib',
        'wheel',
        'torch',
        'numexpr',
        'bottleneck',
        'cuda-python==11.8.1'
    ],
    extras_require={
        'dev': ['h5py'],
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
                'lama_img_info=lama.utilities.lama_img_info:main',
                'lama_ark_imp_pro=lama.scripts.lama_ark_img_pro:main',
                'lama_radiomics_runner=lama.scripts.lama_radiomics_runner:main',
                'lama_two_way_plotter=lama.scripts.two_way_plotter:main',
                'lama_machine_learning=lama.scripts.lama_machine_learning:main'
            ]
        },
)
