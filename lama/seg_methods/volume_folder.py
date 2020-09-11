import SimpleITK as sitk
import argparse
import os

parser = argparse.ArgumentParser(description='get voxel count uCT image')
parser.add_argument('-i', dest='in_dir', help='folder with images in', required=True)
args = parser.parse_args()

vox_file = open(os.path.join(args.in_dir,"voxels.txt"), 'w+')


for line in os.listdir(args.in_dir):
    if line == "voxels.txt":
        continue

    print line
    vox_file.write(line)
    file = os.path.join(args.in_dir,line)

    print "reading file in"
    seg = sitk.ReadImage(file)

    print "calculating volume"
    stats = sitk.StatisticsImageFilter()
    stats.Execute(seg)

    print "Total voxels:"
    vox = stats.GetSum()
    print vox

    vox_file.write("Total voxels: \n")
    vox_file.write(str(vox))

    print "\n"


vox_file.close()


