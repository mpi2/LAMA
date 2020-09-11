import SimpleITK as sitk
from seg import fill_image_all
from seg import fill_image
import sys

img = sitk.ReadImage(sys.argv[1])

print "####First set of hole filling####"
seg = fill_image_all(img)
print "####Second set of hole filling####"
seg = fill_image_all(seg)

sitk.WriteImage(seg,os.path.join(sys.argv[2],"hole_filled.nrrd"))