"""
This module does some post-processing of the stats results
 and writes out the results to file
"""

from pathlib import Path
from typing import Tuple
import numpy as np
from lama.common import write_array

from lama.stats.standard_stats.stats_objects import Stats

MINMAX_TSCORE = 50
FDR_CUTOFF = 0.05


def write(results: Stats,
          mask: np.ndarray,
          root_out_dir: Path,
          stats_name: str):
    """
    Write the line and specimen-level results.

    Remove any Nans
    Threshold the t-statstistics based on q-value
    Write nrrds to file.

    Parameters
    ----------
    results
    out_dir
        The root directory to put the results in
    name
        An the stats type
    """

    # Line-level results
    line = results.input_.line
    line_tstats = results.line_tstats
    line_qvals = results.line_qvals
    shape = results.input_.shape

    line_filt_tstats = result_cutoff_filter(line_tstats, line_qvals)

    line_result = rebuild_array(line_filt_tstats, shape, mask)

    out_dir = root_out_dir / line/ stats_name
    out_dir.mkdir(parents=True, exist_ok=True)

    line_out_path = out_dir / f'{line}_{stats_name}.nrrd'

    write_array(line_result, line_out_path)

    # specimen-level results
    for spec_id , spec_res in results.specimen_results.items():

        spec_t = spec_res['t']
        spec_q = spec_res['q']

        spec_filt_tstats = result_cutoff_filter(spec_t, spec_q)

        spec_result = rebuild_array(spec_filt_tstats, shape, mask)

        spec_out_path = out_dir / f'{spec_id}_{stats_name}.nrrd'

        write_array(spec_result, spec_out_path)


def rebuild_array(array: np.ndarray, shape: Tuple, mask: np.ndarray) -> np.ndarray:
    """
    The stats pipeline uses masked data throughout to save on resources
    This function rebuilds the output files to the orginal 3D sizes of the input volumes

    Parameters
    ----------
    array
        1d masked array to rebuild
    shape
        shape of input volume
    mask
        3D mask

    Returns
    -------
    3d rebuilt array

    """

    array[array > MINMAX_TSCORE] = MINMAX_TSCORE
    array[array < -MINMAX_TSCORE] = - MINMAX_TSCORE

    full_output = np.zeros(shape)
    full_output[mask != False] = array
    return full_output.reshape(shape)


def result_cutoff_filter(t: np.ndarray, q: np.ndarray) -> np.ndarray:
    """
    Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05

    Parameters
    ----------

    """
    if len(t) != len(q):
        raise ValueError
    else:
        mask = q > FDR_CUTOFF
        t[mask] = 0

    return t
