import SimpleITK as sitk

class CastResampleIntensity():
    def __init__(self, pixel_type, resample=False):
        self.pixe_type = pixel_type
        self.resample = resample

    def run(self, inpath, outpath):
        """
        :param inpath:
        :param outpath:
        :raises: IOError
        :return:
        """
        try:
            vol = sitk.ReadImage(inpath)
        except RuntimeError as e:
            raise IOError("Can't open the file {}".format(inpath))
        if self.resample:

            cast = sitk.Cast(sitk.RescaleIntensity(vol), self.pixe_type)
        else:
            cast = sitk.Cast(vol, self.pixe_type)

        sitk.WriteImage(cast, outpath)


class Average():
    def __init__(self):
        pass

