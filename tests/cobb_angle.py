#!/usr/bin/env python

"""

"""



def cobb_angle(infile, output):
    points = parse_vtk(infile)
    thor1, thor2, lumb1, lumb2 = points



def parse_vtk(infile):
    points = []
    with open(infile, 'r') as fh:
        for line in fh:
            try:
                float(line.split()[0])
            except ValueError:
                continue
            x, y, z = [float(x) for x in line.split()]
            points.append((x, z))
    return points

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("icalculate Cobb angles from VTK files")
    parser.add_argument('-i', '--input', dest='input', help='input_vtk', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', help='output dir', required=True)
    args = parser.parse_args()
    cobb_angle(args.input, args.outdir)



    """harry's code

def line3d(x, y, range_=(0, 100)):
    def A(x):
        return np.vstack([x, np.ones(len(x))]).T

    def coeffs(x, y):
        m, c = np.linalg.lstsq(A(x), y)[0]  # defines coefficients y(x) = m * x + c
        return m, c

    c = coeffs(x, y)

    y_ = []
    for x_ in range(*range_):
        y__ = c[0] * x_ + c[1]
        y_.append(y__)
    return np.asarray(y_)

def Cobb_lines(points, image_array=None, *args, **kwargs):

    x, y, z = points.T
    x_range = (0, image_array.shape[2])
    y_range = (0, image_array.shape[1])
    z_range = (0, image_array.shape[0])

    z_ = line3d(x, y, x_range)
    print("z_ = ", z_)
    print("z  = ", z)
    y_ = line3d(z, x, z_range)
    print("y_ = ", y_)
    print("y  = ", y)
    x_ = line3d(z, y, z_range)
    print("x_ = ", x_)
    print("x  = ", x)

    output_array = np.hstack((z_.T, y_.T, x_.T))
    return output_array

def Cobb_angle(points,
               quiet=False, *args, **kwarg):
    _, _, C, D, E, F = points
    CD = C - D
    EF = F - E
    angle = degrees(com.angle_between(CD, EF))
    return angle

def run(array, quiet=False):
    points = get_coordinates(array)
    Cobb_array = Cobb_lines(points)
    angle = Cobb_angle(points)
    print("Cobb angle: {} degrees.\n".format(angle))
    return angle, Cobb_array    # TODO: Rename Cobb_array to something that describes the vars content

def file_filter(root, file, valid_extensions=(".fcsv", ".vtk", ".acsv"), quiet=False):
    file_path = os.path.join(root, file)
    if file.lower().endswith(valid_extensions):  # if input path is a valid file
        if file.lower().endswith(".fcsv"):
            array = com.read_fcsv(file_path, quiet=quiet)
            return run(array, quiet=quiet)
        elif file.lower().endswith(".vtk"):
            array = com.read_vtk(file_path, quiet=quiet)
            return run(array, quiet=quiet)
        else:
            sys.stderr.write("WARNING: File read-in failed. \n"
                             "\t'{1}' is not in a valid format. \n"
                             "\tValid formats: {2}.\n".format(file_path, valid_extensions))
            return
    else:
        # sys.stderr.write("WARNING: File read-in failed. Cannot read '{}/{}'.\n".format(os.path.split(root)[1], file))
        return


def get_Cobb_angles(input_path, quiet=False, *args, **kwargs):
    """ Get Cobb angle. """
    angle_list = []
    if os.path.isfile(input_path):  # if input path is a file
        root, file = os.path.split(input_path)
        angle = file_filter(root, file, quiet=quiet)
        angle_list.append(angle)
    elif os.path.isdir(input_path): # else if input path is a directory
        for root, dirs, files in os.walk(input_path):
            for file in files:
            # try:
                angle = file_filter(root, file, quiet=quiet)
                angle_list.append(angle)
            # except (OSError, IOError, ValueError):  # skip files that aren't in the correct format
            #     sys.stderr.write("WARNING: File '{}' is not in the appropriate format.\n".format(file))
    else:
        sys.stderr.write("ERROR: Input path {} is neither file or directory.\n".format(input_path))
        sys.exit(1)
    angle_list = [angle for angle in angle_list if angle is not None]
    angle_array = np.asarray(angle_list)
    output_path = os.path.join(input_path, "output.csv")
    com.write_csv(output_path, angle_array)  # TODO: Write to file.

if __name__ == '__main__':

    quiet = False

    import argparse

    parser = argparse.ArgumentParser("$")
    parser.add_argument('-i', '--inpath', dest='pathin', help='Input path', required=True)
    parser.add_argument('-o', '--outpath', dest='pathout', help='Output path.', required=True)
    args = parser.parse_args()

    get_Cobb_angles(args.pathin, quiet=quiet)
    """