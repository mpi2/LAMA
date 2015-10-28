"""
Henrik: Reimplementing this to be a hill climber
"""

import numpy as np
cimport numpy as np
import SimpleITK as sitk
import os
import copy


class HcShrink(object):
    def __init__(self, label, label_num, jac_value, padding, out_dir, ngen):
        l = sitk.ReadImage(label)

        cdef float jac_float = float(jac_value)
        bounding_box, array_roi, orig_size = self.get_label_bounding_box(l, label_num, padding, jac_value)

        cdef bb = bounding_box
        self.bb = bb

        cdef np.ndarray ideal = array_roi.astype(np.float32)
        self.ideal = ideal

        # Henrik: not sure what this is for
        cdef list def_shape = list(self.ideal.shape) + [3]
        self.def_shape = def_shape

        self.full_size_def_shape = list(orig_size) + [3]

        self.run(out_dir, int(ngen))

    @staticmethod
    def get_label_bounding_box(label_map, label_num, padding, jac_value):
        shapeXYZ = label_map.GetSize()
        filter = sitk.LabelStatisticsImageFilter()
        filter.Execute(label_map, label_map)
        bb = list(filter.GetBoundingBox(int(label_num))) # in x,y,z

        #add padding
        bb[0] -= padding
        if bb[0] < 0:
            bb[0] = 0
        bb[1] += padding
        if bb[1] > shapeXYZ[0]:
            bb[1] = shapeXYZ[0]

        bb[2] -= padding
        if bb[0] < 0:
            bb[0] = 0
        bb[3] += padding
        if bb[3] > shapeXYZ[1]:
            bb[3] = shapeXYZ[1]

        bb[4] -= padding
        if bb[4] < 0:
            bb[4] = 0
        bb[5] += padding
        if bb[5] > shapeXYZ[2]:
            bb[5] = shapeXYZ[2]

        # Convert bounding box to numpy zyx coords
        bb_numpy = [bb[4], bb[5], bb[2], bb[3], bb[0], bb[1]]


        # test the bounding box
        arr = sitk.GetArrayFromImage(label_map)
        orig_size = arr.shape
        print 'jac', jac_value
        roi_array = arr[bb_numpy[0]:bb_numpy[1], bb_numpy[2]: bb_numpy[3], bb_numpy[4]: bb_numpy[5]].astype(np.float32)
        print roi_array.max(), roi_array.min()
        roi_array[roi_array == label_num] = jac_value
        print roi_array.max(), roi_array.min()
        roi_array[roi_array != jac_value] = 1.0
        print roi_array.min(), roi_array.max()

        return bb_numpy, roi_array, orig_size




    @staticmethod
    # get the determinents per voxel
    def make_jac(np.ndarray individual):

        im = sitk.GetImageFromArray(individual)
        jac = sitk.DisplacementFieldJacobianDeterminant(im)
        cdef np.ndarray jac_array = sitk.GetArrayFromImage(jac)
        return jac_array

    # calculate the fitness
    def calc_fitness(self, np.ndarray individual):
        """
        Add weight to regions where ideal_jac != 1.0
        """
        jac_array = self.make_jac(individual)
        cdef float surronding_comp = np.sum(np.square(self.ideal[self.ideal == 1.0] - jac_array[self.ideal == 1.0])) / individual.size
        cdef float organ_comp = np.sum(np.square(self.ideal[self.ideal != 1.0] - jac_array[self.ideal != 1.0])) / individual.size

        return surronding_comp + (organ_comp *3)

    # evalualte a whole population
    def evaluate(self, np.ndarray ind):
        return self.calc_fitness(ind)

    def mutate(self, np.ndarray ind, float fit):

        # pick a random position
        cdef int random_position = np.random.randint(0, ind.size)
        cdef float mutate_val = np.random.uniform(0, fit *  30)
        cdef float rand = np.random.uniform(-mutate_val, mutate_val)
        ind.ravel()[random_position] += rand
        return ind

    def vector_to_def(self, np.ndarray vector):
        cdef np.ndarray reshaped_ind  = np.array(vector).reshape(self.def_shape)
        return reshaped_ind

    def get_initial_individual(self, shape, val):

        cdef np.ndarray ind = np.zeros(shape, dtype=np.float32)
        #ind = np.random.uniform(-val, val, shape).astype(np.float32)
        #ind[self.ideal == 1.0] = 0.0
        return ind

    def run(self, out_dir, int ngen):

        #Print the ideal jac
        ideal_out = os.path.join(out_dir, 'ideal_jac.nrrd')
        ideal_im = sitk.GetImageFromArray(self.ideal)
        sitk.WriteImage(ideal_im, ideal_out)
        cdef tuple ideal_shape = self.ideal.shape

        cdef list vector_field_shape = list(ideal_shape)
        vector_field_shape.append(3)  # 2 for 2d vectors

        print vector_field_shape

        cdef np.ndarray ind = self.get_initial_individual(vector_field_shape, 0.0)

        #print ind.shape

        temp_results = os.path.join(out_dir, 'intermediate_results')
        if not os.path.exists(temp_results):
            os.mkdir(temp_results)

        cdef int intermediate_num = 3000  # Spit out intermediate results every x num generations
        cdef float fit

        cdef bint temp
        cdef np.ndarray  temp_ind

        with open(os.path.join(out_dir, 'chart_data.txt'), 'w')as dh:

            for gen in range(ngen):

                # evaluate individual
                fit = self.evaluate(ind)

                #print fit

                temp = True
                # mutate the individual until it gets better
                while temp:
                    temp_ind = copy.deepcopy(ind)
                    temp_ind = self.mutate(temp_ind, fit)

                    if self.evaluate(temp_ind) < fit:
                        ind = temp_ind
                        temp = False
                dh.write('{}\n'.format(fit))
                dh.flush()

                if gen != 0 and gen % intermediate_num == 0:
                    self.write_results(temp_results, ind)
            self.write_results(out_dir, ind)

    def write_results(self, out, ind):

        jac_result = self.make_jac(ind)

        # Insert result into

        # write the jac
        outpath = os.path.join(out, "jac_{}.nrrd".format(1))
        sitk.WriteImage(sitk.GetImageFromArray(jac_result), outpath, True)

        # write the def
        outpath = os.path.join(out, "def_{}.nrrd".format(1))
        def_roi = self.vector_to_def(ind)
        big_def = np.zeros(np.prod(self.full_size_def_shape)).reshape(self.full_size_def_shape)
        big_def[self.bb[0]: self.bb[1], self.bb[2]: self.bb[3], self.bb[4]: self.bb[5]] = def_roi


        sitk.WriteImage(sitk.GetImageFromArray(big_def), outpath, True)





if __name__ == '__main__':
    import sys

    label_map = sys.argv[1]
    label_num = sys.argv[2]
    jac_value = sys.argv[3]
    padding = sys.argv[4]
    out_dir = sys.argv[5]
    ngen = sys.argv[6]

    HcShrink(label_map, label_num, jac_value, padding, out_dir, ngen)


