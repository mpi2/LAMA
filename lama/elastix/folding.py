"""
Folding in the deformations cn casue problems with the inversions.
This module corrects overafolding and ensures injectivity of the transform
"""
import numpy as np 
from pathlib import Path
from typing import Union


K2 = 2.046392675
K3 = 2.479472335
A2 = np.sqrt(
    (3/2)**2 + (K2 - (3/2))**2
)
A3 = np.sqrt(
    ((3/2)**2 + (K2 - (3/2))**2 + (K3 - K2)**2)
)


class BSplineParse():

    def __init__(self, tform_file: str):
        """
        Read in the coefficients from an elastix tform parameter files

        Parameters
        ----------
        tform_file
            the path to the elastix TranformParameter file

        Attributes
        ----------
        coefs: np.ndarray
            m*n array. m number of control points, n = num axes

        Returns
        -------

        tuple

        0:
            if Bspline:
                n*m np.array where n is number of control points and m is vector of coordinates (xyz)
            if Euler, similarity, affine
                1d array of transormation parameters
        1:
            list of non-transform parameter lines from the tform file

        """

        with open(tform_file, 'r') as fh:

            other_data = []

            bspline = False

            for line in fh:
                if line.startswith('(TransformParameters '):
                    tform_str = line.strip()
                    tform_str = tform_str.strip(')')
                    tform_str = tform_str.lstrip('(TransformParameters ')
                    tform_params = tform_str.split(' ')

                    tform_params= [float(x) for x in tform_params]
                else:
                    other_data.append(line)

                if line.startswith('(Transform "BSplineTransform")'):
                    bspline = True

            if bspline:
                tform_params = np.array(np.array_split(tform_params, 3)).T

            self.coefs = tform_params
            self.elastix_params = other_data

    def control_point_coords(self):

        for line in self.elastix_params:

            if line.startswith('(Size '):
                size = [int(x) for x in line.strip().split(' ')[1:]]

            elif line.startswith('GridOrigin'):
                grid_origin = [int(x) for x in line.strip().split(' ')[1:]]

            elif line.startswith('GridSpacing'):
                grid_spacing = [int(x) for x in line.strip().split(' ')[1:]]


def condition_1(coefs):
    """
    numpy array of tform coefs of shape n_coefs * dims

    Returns
    -------

    """
    t = np.abs(coefs) < 1/K3
    condition_met = np.apply_along_axis(np.all, 1, t)
    return condition_met


def condition_2(coefs):
    # def f(x):
    #     # s = np.sum(x) < (1 / A3) ** 2
    #     s = np.linalg.norm(x) < (1 / A3)
    #     return s
    # t = coefs**2
    # condition_met = np.apply_along_axis(f, 1, t)
    # return condition_met
    return np.linalg.norm(coefs, axis=1) < (1 / A3)


# def correct_1(coefs, potential_fold_indices):
#     for idx in potential_fold_indices:
#         bounds = (1 / A3) ** 2
#         coefs[idx, :] = bounds
#     return coefs


def correct(coefs):
    coefs = np.copy(coefs)

    result = []
    for c in coefs:

        # Make this vector satisfy condition 2
        unfolded = np.array(c) / np.linalg.norm(c) * (1 / A3)
        # unfolded = c * (max_bound / np.sum(c))

        # Replace the vector of coefs
        result.append(unfolded)

    return result


def unfold_bsplines(bs: Union[BSplineParse, str], outfile=None) -> np.ndarray:

    if isinstance(bs, (str, Path)):
        bs = BSplineParse(bs)

    # Scale the coeficents by the grid size
    spacing_str = bs.elastix_params[21].split(' ')[1:4]
    grid_size = [float(x.strip().strip(')')) for x in spacing_str]

    # this conversion may be needed if bs.coefs doesn't accept numpy array operations
    # scale to (1, 1)-spacing grid

    bs.coefs[:, 0] /= grid_size[0]
    bs.coefs[:, 1] /= grid_size[1]
    bs.coefs[:, 2] /= grid_size[2]
    # grid points that need to be corrected
    idx_to_correct = ~(condition_1(bs.coefs) | condition_2(bs.coefs))

    if not np.any(idx_to_correct): # No folding so write or return orginal
        if outfile:
            write_tform(bs.coefs, bs.elastix_params, outfile)
            return
        else:
            return bs
    # correct potentially problematic grid points
    bs.coefs[idx_to_correct, :] = correct(bs.coefs[idx_to_correct, :])
    # rescale to original grid spacing
    bs.coefs[:, 0] *= grid_size[0]
    bs.coefs[:, 1] *= grid_size[1]
    bs.coefs[:, 2] *= grid_size[2]

    if outfile:
        write_tform(bs.coefs, bs.elastix_params, outfile)

    return bs.coefs


def write_tform(tform_params, other_data, tform_file, init_tform=None):
    """
    Write elastix transform file

    Parameters
    ----------
    tform_params
        The spline coefficient or tranform parameters
        m*n array. m = num control points, n = dims of coefs (3)
    other_data
        Lines of the other data from the transform parameer file

    tform_file
        Where to write the tform

    init_tform

    Returns
    -------

    """

    with open(tform_file, 'w') as fh:

        for line in other_data:

            if init_tform and line.startswith('(InitialTransformParametersFileName'):
                line = f'(InitialTransformParametersFileName "{init_tform}")\n'

            if line.startswith('(NumberOfParameters '):
                num_params = int(line.strip().strip(')').split(' ')[1])

            fh.write(line)

        if tform_params.size != num_params:
            raise ValueError(f"stated num params: {num_params} does not match actual {tform_params.size}")
        # write out all the x column first by using fortran mode. Set print options to prevent ... shortening
        np.set_printoptions(threshold=np.prod(tform_params.shape), suppress=True)

        fh.write(f"(TransformParameters ")
        tform_params.flatten(order='f').tofile(fh, sep=" ", format="%.6f")
        fh.write(')')
