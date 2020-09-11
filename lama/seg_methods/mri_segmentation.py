import SimpleITK as sitk
import os, sys
import subprocess
import shutil
import argparse


def get_init_seed(img,output):
    print "########## First stage: Getting the Seed ##########"
    print "threshold"
    output = os.path.join(output,"initial_seed")

    if not os.path.exists(output):
        os.makedirs(output)

    seg = sitk.BinaryThreshold(img, lowerThreshold = 900, upperThreshold = 20000, insideValue = 255, outsideValue = 0 )
    sitk.WriteImage(seg,os.path.join(output,"threshold.nrrd"))

    print "get the correct component"
    seg = sitk.ConnectedComponent(seg)
    sitk.WriteImage(seg,os.path.join(output,"componet_all.nrrd"))
    # change the order of labels to largest first
    seg = sitk.RelabelComponent(seg)
    #seg = seg == 5202
    seg = seg == 1
    sitk.WriteImage(seg,os.path.join(output,"component_largest.nrrd"))

    print "erode"
    for x in range(0, 8):
        print "Erode number %d" % (x)
        seg =  sitk.BinaryErode(seg != 0)


    print "writing erode"
    sitk.WriteImage(seg,os.path.join(output,"erode.nrrd"))

    print "getting component"
    seg = sitk.ConnectedComponent(seg)
    seg = sitk.RelabelComponent(seg)
    sitk.WriteImage(seg,os.path.join(output,"component_eroded.nrrd"))
    seg = seg == 1

    sitk.WriteImage(seg,output+"seed1.nrrd")
    print "Dilate the seg"
    for x in range(0, 3):
        print "Dilate number %d" % (x)
        seg =  sitk.GrayscaleDilate(seg,10)
    sitk.WriteImage(seg,os.path.join(output,"seed_dilated.nrrd"))
    return seg


def level_set_seg(img,seed,output):
    print "######## Second stage: Perform level set threshold segmentation ###########"
    output = os.path.join(output,"level_set_seg")

    if not os.path.exists(output):
        os.makedirs(output)

    print "get stats of image"
    #stats = sitk.LabelStatisticsImageFilter()
    #stats.Execute(img, seg)

    print "get threshold limits"
    factor = 3
    lower_threshold = 200 #stats.GetMean(1)-factor*stats.GetSigma(1)
    upper_threshold = 20000 #stats.GetMean(1)+factor*stats.GetSigma(1)

    print "get distance map"
    seed_map = sitk.SignedMaurerDistanceMap(seed, insideIsPositive=True, useImageSpacing=True)
    sitk.WriteImage(seed_map,os.path.join(output,"initial_seed.nrrd"))

    print "set up threshold parameters"
    lsFilter = sitk.ThresholdSegmentationLevelSetImageFilter()
    lsFilter.SetLowerThreshold(lower_threshold)
    lsFilter.SetUpperThreshold(upper_threshold)
    lsFilter.SetMaximumRMSError(0.02)
    lsFilter.SetNumberOfIterations(10000)
    lsFilter.SetCurvatureScaling(1)
    lsFilter.SetPropagationScaling(1)
    #lsFilter.ReverseExpansionDirectionOn()

    print "Paramaters:"
    print lsFilter

    print "change original image to float32"
    img_f = sitk.Cast(img, sitk.sitkFloat32)
    sitk.WriteImage(img_f,os.path.join(output,"img_f.nrrd"))

    niter = 0
    for i in range(0, 3):
        print "perform ls seg"
        ls = lsFilter.Execute(seed_map, img_f)

        sitk.WriteImage(ls,os.path.join(output,str(i)+"_ls.nrrd"))
        niter += lsFilter.GetNumberOfIterations()
        print "LevelSet after "+str(niter)+" iterations and RMS "+str(lsFilter.GetRMSChange())
        print "get connected components"
        connected = sitk.ScalarConnectedComponent(ls)

        sitk.WriteImage(connected,os.path.join(output,"ls_compon"+str(i)+".nrrd"))

        seg = connected > 1

        sitk.WriteImage(seg,os.path.join(output,"ls_seg"+str(i)+".nrrd"))

        seed_map = sitk.SignedMaurerDistanceMap(seg, insideIsPositive=True, useImageSpacing=True)

    sitk.WriteImage(connected,os.path.join(output,"ls_compon.nrrd"))
    sitk.WriteImage(seg,os.path.join(output,"ls_seg.nrrd"))

    return seg

def clean_up(img,seg,output):
    print "######## Third stage: Clean up ###########"
    output = os.path.join(output,"clean_up")

    if not os.path.exists(output):
        os.makedirs(output)

#    print "dilate"
#    seg = sitk.BinaryDilate(seg,30)
#     dilate = sitk.DilateObjectMorphologyImageFilter()
#     dilate.SetKernelRadius(5)
#     dilate.execute(seg)

    sitk.WriteImage(seg,os.path.join(output,"dilate.nrrd"))

    print "fill hole"
    seg = sitk.BinaryFillhole(seg)
    sitk.WriteImage(seg,os.path.join(output,"ls_seg_hole_fill.nrrd"))

    print "mask"
    masked_image1 = sitk.Mask(img, seg)
    sitk.WriteImage(masked_image1,os.path.join(output,"masked1.nrrd"))

    seg = sitk.BinaryThreshold(masked_image1, lowerThreshold = 0, upperThreshold = 20000, insideValue = 255, outsideValue = 0 )
    masked_image2 = sitk.Mask(masked_image1, seg)
    sitk.WriteImage(masked_image2,os.path.join(output,"masked2.nrrd"))

    #seg = sitk.GradientMagnitudeRecursiveGaussian(masked_image2, sigma=0.01)
    #sitk.WriteImage(seg,os.path.join(output,"blur.nrrd"))
    #sitk.WriteImage(seg,os.path.join(outputFolder,filter+"_blur.nrrd"))
    # Perform the watershed alg
    print "Doing watershed"
    seg = sitk.MorphologicalWatershed(masked_image2, level=0.9, markWatershedLine=False, fullyConnected=False)
    sitk.WriteImage(seg,os.path.join(output,"ws.nrrd"))
    # This makes the muliple componets of the watershed segmentation one
    print "connecting components"
    seg = sitk.ConnectedComponent(seg != seg[0, 0, 0])
    seg = sitk.RelabelComponent(seg)
    seg = seg == 1
    sitk.WriteImage(seg,os.path.join(output,"components.nrrd"))

    # The watershed on its own doesnt quite get everything. If we dilate it catches the missing bits. Could increase the dilation if want to be sure
    # random bits are caught
    # Has to be greater than 1
    masked_image3 = sitk.Mask(masked_image2, seg)
    sitk.WriteImage(masked_image3,os.path.join(output,"masked3.nrrd"))

    print "filling holes"
    seg = sitk.BinaryFillhole(seg!=0)

    sitk.WriteImage(masked_image3,os.path.join(output,"seg3_embryo.nrrd"))

    print "overlay"
    tif = os.path.join(output,"overlay_ls.tif")
    if os.path.exists(tif):
        os.remove(tif)
    overlay = sitk.LabelOverlay(img, seg)
    sitk.WriteImage(overlay,os.path.join(output,"overlay_ls.tif"))





parser = argparse.ArgumentParser(description='Method of extracting embryos from the MRI dataset')
parser.add_argument('-i', dest='in_dir', help='mri image', required=True)
parser.add_argument('-o', dest='out_dir', help='destination for output', required=True)
parser.add_argument('-s', dest='seed', help='initial seed if already known', required=False)
parser.add_argument('-d', dest='dilate', help='dilate if different from default', required=False)
parser.add_argument('-e', dest='erode', help='erode if different from default', required=False)
parser.add_argument('-l', dest='level_set', help='level set seg if different from default', required=False)
args = parser.parse_args()

input = args.in_dir
output = args.out_dir

if not os.path.exists(output):
    os.makedirs(output)

print "reading image"
img = sitk.ReadImage(input)

# Get initial seed
if args.level_set:
    print "no seed required"
elif args.seed:
    seed = sitk.ReadImage(args.seed)
else :
    seed = get_init_seed(img,output)

# perform level set
if args.level_set:
    print "using level set seg provided"
    seg = sitk.ReadImage(args.level_set)
else:
    seg = level_set_seg(img,seed,output)

# perform clean up
clean_up(img,seg,output)





