from pathlib import Path
from typing import Tuple


def final_red_cyan_iterator(root, orientation='sagittal') -> Tuple[Path, str, str]:
    """
    Get the red/cyan overlay of the final registered image

    Parameters
    ----------
    root
        Directory containing one or more lama lines directories
    orientation
        sagittal, coronal or axial
    Returns
    -------
    0: Path to the image
    1: line_id
    2: specimen id
    """
    for line_dir in root.it