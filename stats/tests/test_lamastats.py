import nose
from nose import with_setup
from staging.get_volumes_by_stage import VolumeGetter
from tempfile import NamedTemporaryFile

"""
nose tests for bu_lama_stats.py
bu_lama_stats.py is the main script for running the statistical analysis for LAMA
"""



wt_staging_file = NamedTemporaryFile()
mut_staging_file = NamedTemporaryFile()


def setup():
    print 'setting up'
    with open(wt_staging_file.name, 'w') as fh:
        fh.write(wt_data)

def save_mutant_file(data):
    with open(mut_staging_file.name, 'w') as mf:
        mf.write(data)


@with_setup(setup)
def test_get_config():
    """
    Test the function that gets all the
    Returns
    -------

    """
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    #stager.plot() # Coud write plot file
