

class Registration(object):
    def __init__(self):
        pass

    def set_parameter_file(self, param):
        """
        For elastix, this will be the parameter file to e read in by elastix
        For tools that are cmd line driven, the parameters will be read and fed to cmd line
        """

    def set_fixed_image(self, img):
        raise NotImplementedError

    def set_moving_image(self, img):
        raise NotImplementedError

class ElastixRegistration(Registration):
    