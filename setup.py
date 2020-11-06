# coding: utf-8
import platform
import subprocess
import os
import re
import sys

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake','--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))
        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                                   out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required for Windows")
        for ext in self.extensions:
            self.build_extension(ext)


def build_extension(self, ext):
    extdir = os.path.abspath(
        os.path.dirmane(self.get_ext_fullpath(ext.name)))
    cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir, '-DPYTHON_EXECUTABLE=' + sys.executable]

    cfg = 'Debug' if self.debug else 'Release'
    build_args = ['--config', cfg]

    if platform.system() == "Windows":
        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
            cfg.upper(),
            extdir)]
        if sys.maxsize > 2**32:
            cmake_args += ['-A', 'x64']
        build_args += ['--', '/m']
    else:
        cmake_args += ['-DCMAKE_BUILD_TYPE' + cfg]
        build_args += ['--','-j2']

    env = os.environ.copy()
    env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\'.format(
        env.get('CXXFLAGS' ''),
        self.distribution.get_version())
    if not os.path.exists(self.build_temp):
        os.makedirs(self.build_temp)
    subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                          cwd=self.build_temp, env=env)
    subprocess.check_call(['cmake','--build', '.'] + build_args,
                          cwd=self.build_temp)
    print()



setup(
    name='lama_phenotype_detection',
    download_url='https://github.com/mpi2/lama/archive/0.9.4.tar.gz',
    version='0.9.53',
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

setup(
    name='mesher',
    version='0.1',
    author='Kyle_Drover',
    author_email='kyle.drover@anu.edu.anu.au',
    description='imports ITKSNAP cluster and snake contouring',
    long_description='',
    packages=find_packages('src'),
    package_dir={'':'src'},
    ext_modules=[CMakeExtension('mesher/itksnap_mesher')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
