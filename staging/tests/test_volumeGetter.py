import nose
from nose import with_setup
from staging.get_volumes_by_stage import VolumeGetter
from tempfile import NamedTemporaryFile

# Need to simulate files being passed to stage file getter
# Not sure how to do that yet so I'll use a tempfile object for now
#  wt_staging_file, mut_staging_file, littermate_basenames=None, plot_path=None



wt_staging_file = NamedTemporaryFile()
mut_staging_file = NamedTemporaryFile()

wt_data = """vol,value
a,1
b,2
c,3
d,4
e,5
f,6
g,7
h,8
i,9
j,10
k,11
l,11.5
m,12.0
"""


def setup():
    print 'setting up'
    with open(wt_staging_file.name, 'w') as fh:
        fh.write(wt_data)

def save_mutant_file(data):
    with open(mut_staging_file.name, 'w') as mf:
        mf.write(data)


@with_setup(setup)
def test_get_all_wt_within_mut_range():
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    # stager.plot() # Coud write plot file


@with_setup(setup)
def test_removal_of_littermates():
    """
    test whether littermates file works.
    Littermates csv gives IDs of wild type littermates. These should be discounted from the staging calculation
    as they are often a lot larger and can interfere with selection of wild type
    """
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
littermate1,12.0"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name, littermate_basenames=['littermate1'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']


def test_removal_of_littermates_without_extension():

# See whether we can just use basenames instead of full paths
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
littermate1.nrrd,12.0"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name, littermate_basenames=['littermate1'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']


@with_setup(setup)
def test_out_of_range():
    """
    Try with mutants that are outside of the specified range of wild types.
    If not able to get min number of wts, should return None
    """
    mut_data = """vol,value
mut1,13
mut2,14
mut3,15"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files is None


@with_setup(setup)
def test_out_of_range_with_constraint_removed():
    """
    Try with mutants that are outside of the specified range of wild types, but with no constraint on the rnage of sizes
    """
    mut_data = """vol,value
mut1,13
mut2,14
mut3,15"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids(ignore_constraint=True)
    assert files == ['f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']

    mut_data = """vol,value
mut1,0.1
mut2,0.2
mut3,0.3"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids(ignore_constraint=True)
    assert files == ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


@with_setup(setup)
def test_mutant_ids():
    """
    Stager takes mut_ids argument. This is used when only a subset of the mutants are to be used. If not used all
    the mutants in the staging file are used
    """
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
mut_to_ignore,12"""
    save_mutant_file(mut_data)
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name, mut_ids=['mut1', 'mut2', 'mut3'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

    # If no mut_ids are specified, use all available mutants, which should give more wild types back
    stager = VolumeGetter(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']