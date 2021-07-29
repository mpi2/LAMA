from skimage.measure import regionprops
from pathlib import Path
from dataclasses import dataclass
from typing import Union, List, Iterator
import nrrd
import numpy as np

@dataclass
class Slice:
    axis: int
    indices: np.lib.index_tricks.IndexExpression
    image: np.ndarray = None
    labels: np.ndarray = None


rotations = {
    # Number of rotations for each axis to get into our standard orientation
    0: 1,
    1: 1,
    2: 3
}


class AtlasProps:
    def __init__(self, atlas_path: Union[Path, str],
                 raw_image_path: Union[Path, str] = None):

        self.atlas_path = atlas_path
        self.atlas, self.atlas_head = nrrd.read(atlas_path)
        self.atlas = self.atlas.astype(np.int)

        if raw_image_path:
            self.image_path = raw_image_path
            self.image, self.image_head = nrrd.read(raw_image_path)
        else:
            self.image = self.image_head = None

        # Created a dict of regions props indexed by label and sorted by area
        self.props = {prop['label']: prop for prop in sorted(regionprops(self.atlas), key=lambda x: x.area, reverse=True)}

    def slice_indices(self, label: int, axis=0, pad=0) -> List[np.lib.index_tricks.IndexExpression]:
        # get a list of np.s_ slices which can be used to generate 2D slices form the atlas
        bbox = self.props[label].bbox

        slices = []

        # Test
        start = bbox[2] - pad
        end = bbox[5] + pad

        for i in range(start, end):
            slices.append(np.s_[bbox[0]-pad: bbox[3]+pad, bbox[1]-pad: bbox[4]+pad, i])

        return slices

    # def slices(self, label, n_slices, axis=0, pad=0, return_labels=True, return_image=True)  -> Iterator[slice]:
    #
    #     if not return_image and not return_labels:
    #         raise ValueError('return_labels or return_image must be True')
    #
    #     indices = self.slice_indices(label, axis, pad)
    #
    #     s = Slice()
    #
    #     for idx in indices:
    #         s.indices = idx
    #         s.axis = axis
    #         if return_image:
    #             s.image = self.image[idx]
    #         if return_labels:
    #             s.labels = self.atlas[idx]
    def get_roi_slices(self, label, n_slices, axis, padding, main_axis_padding=3, dtype=None):
        """
        slice_padding
            The ammount of padding applied to each 2D slice
        slice_padding_2
            The ammount of padding used for the indices of the main axis
        """
        for roi in self.get_roi(label, padding):
            m = padding - main_axis_padding
            if m < padding:
                m = padding

            slices = []
            main_axis_len = roi.shape[axis]
            slice_indices = list(np.linspace(0 + m, main_axis_len - m, n_slices, dtype=int))

            for idx in slice_indices:

                s = np.rot90(np.take(roi, idx, axis=axis), k=rotations[axis])
                if dtype:
                    s = s.astype(dtype)
                slices.append(s)
            yield slices

    def get_roi(self, label, padding: int=0, shape=None) -> np.ndarray:

        if not shape: # Use the label bounding box
            b = self.props[label].bbox

        else:  # Use centroid and shape
            c = self.props[label].centroid
            b = [int(x) for x in [c[0] - (shape[0] / 2), c[1] - (shape[1] / 2), c[2] - (shape[2] / 2),  # bbox starts
                 c[0] + (shape[0] / 2), c[1] + (shape[1] / 2), c[2] + (shape[2] / 2)] ] # bbox ends

        if padding:
            b = [x - padding if i < 3 else x + padding for i, x in enumerate(b)]
            b = np.clip(b, 0, max(b))
        atlas_roi = self.atlas[b[0]:b[3], b[1]: b[4], b[2]: b[5]]
        image_roi = self.image[b[0]:b[3], b[1]: b[4], b[2]: b[5]]

        atlas_roi[atlas_roi != label] = 0

        return atlas_roi, image_roi


def remove_unconected(image, label_num) -> np.ndarray:
    """
    Keep only the largest coonected componenet label with given label num. Zero out all other labels

    """
    from skimage.measure import label, regionprops
    image = np.copy(image)
    image[image != label_num] = 0
    relab = label(image)
    biggest_label = sorted(regionprops(relab), key=lambda x: x.area, reverse=True)[0].label
    image[:] = 0
    image[relab == biggest_label] = label_num
    return image