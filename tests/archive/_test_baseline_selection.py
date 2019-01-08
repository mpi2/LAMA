"""
Test the BaelineSelector. The calss that given a bunch of bselines and mutants along with some staging metric such as
crown-rump length (CRL) return a list of baselines to use in the analysis of a given line
"""


from nose import with_setup
from staging.baseline_selection import BaselineSelector
from tempfile import NamedTemporaryFile
from nose.tools import nottest


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
    with open(wt_staging_file.name, 'w') as fh:
        fh.write(wt_data)


def save_mutant_file(data):
    with open(mut_staging_file.name, 'w') as mf:
        mf.write(data)


@with_setup(setup)
def test_exclude_small_mutants():
    """
    Mutants that are too small/big should be removed, and should not contribute to the baseline selection
    Returns
    -------

    """
    mut_data = """vol,value
mut1,0.2
mut2,2.0
mut3,10.0"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, mut_ids=['mut1', 'mut2', 'mut3'])
    exclude_muts = stager.excluded_mutants
    assert (exclude_muts == ['mut1'])

    # This breaks as the mutants are took out in run_lama_stats. Fix this in next iteration
    files = stager.filtered_wt_ids()
    assert files == ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'] # If mut1 was being included we would get baselines down to 'a' also
#
@with_setup(setup)
def test_exclude_large_mutants():
    """
    Mutants that are too small/big should be removed, and should not contribute to the baseline selection
    Returns
    -------

    """
    mut_data = """vol,value
mut1,2.0
mut2,9.0
mut3,15.0"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, mut_ids=['mut1', 'mut2', 'mut3'])
    exclude_muts = stager.excluded_mutants
    assert (exclude_muts == ['mut3'])

    # This breaks as the mutants are took out in run_lama_stats. Fix this in next iteration
    files = stager.filtered_wt_ids()
    assert files == ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']


@with_setup(setup)
def test_get_all_wt_within_mut_range():
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    excluded = stager.excluded_mutants
    # stager.plot() # Coud write plot file


@with_setup(setup)
def test_removal_of_littermates():
    """
    Littermate wild types should be discounted from the baseline set if they are outside the range of the mutants.

    """
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
littermate1,12.0"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, littermate_basenames=['littermate1'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    littermates_to_use = stager.littermates_to_include()
    assert littermates_to_use is None  # Too big


@with_setup(setup)
def test_retain_littermates():
    """
    Littermate wild types should be included in the baseline set if they are outside the range of the mutants.
    """
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
littermate1,9.0"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, littermate_basenames=['littermate1'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

    littermates_to_use = stager.littermates_to_include()
    assert littermates_to_use[0] == 'littermate1'


def test_removal_of_littermates_with_extension():
    """See whether we can just use extensions on ids"""
    mut_data = """vol,value
mut1,3.0
mut2,5.0
mut3,10.0
littermate1.nrrd,12.0"""
    save_mutant_file(mut_data)
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, littermate_basenames=['littermate1'])
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
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert not files


# @with_setup(setup)
# def test_out_of_range_with_constraint_removed():
#     """
#     Try with mutants that are outside of the specified range of wild types, but with no constraint on the rnage of sizes
#     we should get the 8 (this is hard coded in staging) smallest or largest
#     """
#     mut_data = """vol,value
# mut1,13
# mut2,14
# mut3,15"""
#     save_mutant_file(mut_data)
#     stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name)
#     files = stager.filtered_wt_ids(ignore_constraint=True)
#     assert files == ['f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
#
#     mut_data = """vol,value
# mut1,0.1
# mut2,0.2
# mut3,0.3"""
#     save_mutant_file(mut_data)
#     stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name)
#     files = stager.filtered_wt_ids(ignore_constraint=True)
#     assert files == ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


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
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name, mut_ids=['mut1', 'mut2', 'mut3'])
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

    # If no mut_ids are specified, use all available mutants, which should give more wild types back
    stager = BaselineSelector(wt_staging_file.name, mut_staging_file.name)
    files = stager.filtered_wt_ids()
    assert files == ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
#
