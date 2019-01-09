from distutils.core import setup

setup(
    name='lama',
    version='0.9',
    packages=['', 'lama', 'lama.qc', 'lama.dev', 'lama.dev.validation', 'lama.dev.organ_shrinking', 'lama.lib',
              'lama.stats', 'lama.stats.dev', 'lama.stats.old', 'lama.stats.standard_stats',
              'lama.stats.permutation_stats', 'lama.elastix', 'lama.staging', 'lama.img_processing',
              'lama.registration_pipeline', 'tests', 'tests.archive', 'scripts', 'utilities'],
    url='https://www.har.mrc.ac.uk/',
    license='Apache2',
    author='neil',
    author_email='bit@har.mrc.ac.uk',
    description='Phenotype detection pipeline'
)
