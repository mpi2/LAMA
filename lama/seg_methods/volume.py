import SimpleITK as sitk
import argparse
import os

parser = argparse.ArgumentParser(description='get voxel count uCT image')
parser.add_argument('-bin', dest='bin', help='binary segmentated image', required=True)
parser.add_argument('-img', dest='img', help='original image', required=False)
args = parser.parse_args()

seg = sitk.ReadImage(args.bin)
# img = sitk.ReadImage(args.img)

print "calculating volume"
stats = sitk.StatisticsImageFilter()

stats.Execute(seg)

print "Total voxels:"
print stats.GetSum()

# print "get stats "
# stats = sitk.LabelStatisticsImageFilter()
# stats.Execute(img,seg)
# print stats

# print "get shape stats "
# shape_stats = sitk.LabelShapeStatisticsImageFilter()
# shape_stats.Execute(seg)
# print shape_stats


