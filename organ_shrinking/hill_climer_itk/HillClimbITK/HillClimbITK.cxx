/* 
 * File:   main.cpp
 * Author: neil
 *
 * Created on 17 December 2015, 12:00
 */
#include "itkImage.h"
#include "itkRandomImageSource.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream  
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>

using namespace std;

typedef itk::Image<float, 3> FloatImageType;
typedef itk::Image<char, 3> ByteImageType;
typedef itk::Image<float, 4> VectorImageType;

/*
 * 
 */


class Shrinker
{
    
public:
    ByteImageType::SizeType boundingBox;
    FloatImageType::Pointer idealJacImage; 
    VectorImageType::SizeType defImageShape;
    VectorImageType::SizeType fullSizeDefShape;
    Shrinker(string labelPath, int labelNum, float jacValue, int padding, string outDir, int numGens);
    void getLabelBoundingBox(ByteImageType::Pointer labelMap, int labelNumber, int padding, float jacobianvalue);
    
};


Shrinker::Shrinker(string labelPath, int labelNum, float jacValue, int padding, string outDir, int numGens){
    
    typedef itk::ImageFileReader<ByteImageType> ByteReaderType;
    ByteReaderType::Pointer reader=ByteReaderType::New();
    reader->SetFileName(labelPath);
    reader->Update();
    ByteImageType::Pointer labelImage=reader->GetOutput();
    
    ByteImageType::SizeType originalSize;
    originalSize = labelImage->GetLargestPossibleRegion().GetSize();
    
    fullSizeDefShape[0] = originalSize[0];
    fullSizeDefShape[1] = originalSize[1];
    fullSizeDefShape[2] = originalSize[2];
    fullSizeDefShape[3] = 3;
    
    std::cout << "full size def" << fullSizeDefShape << std::endl;
    getLabelBoundingBox(labelImage, labelNum, padding, jacValue);
 }



void Shrinker::getLabelBoundingBox(ByteImageType::Pointer labelMap, int labelNumber, int padding, float jacobianvalue ){
    // Set the padded bounding box, make the ideal jacobian, 
    
    typedef itk::LabelStatisticsImageFilter<typename ByteImageType::Pointer,typename ByteImageType> LabelFilterType;
    LabelFilterType::Pointer labelFilter = new LabelFilterType::New();
    labelFilter->
    
}



int main(int argc, char** argv) {

    if(argc < 6)
    {
    std::cerr << "Usage:\n" << "label map path<string>\nlabel number<int>"
            "\ndesired jacobian<float>\npadding<int>\nout folder<string>\nnumber of generations<int>" << std::endl;
    return EXIT_FAILURE;
    }
    
    std::string labelMapPath = argv[1];
    
    //int labelNumber = atoi(argv[2]);
    std::istringstream iss2( argv[2] );
        int labelNum;

        if (iss2 >> labelNum)
        {
            std::cout << labelNum << std::endl;
        }else
        {
            std::cout << "Need to pass a number for label" << std::endl;
            return 1;
        }
    
    std::istringstream iss3( argv[3] );
        float jacobianValue;

        if (iss3 >> jacobianValue)
        {
            std::cout << jacobianValue << std::endl;
        }else
        {
            std::cout << "Need a float value for the desired jacobian" << std::endl;
            return 1;
        }
    
    std::istringstream iss4( argv[4] );
        int padding;

        if (iss4 >> padding)
        {
            std::cout << padding << std::endl;
        }else
        {
            std::cout << "Need an int value for the padding amount" << std::endl;
            return 1;
        }
    
    std::string outDir = argv[5];
     
    std::istringstream iss6( argv[6] );
        int numGens;

        if (iss6 >> numGens)
        {
            std::cout << numGens << std::endl;
        }else
        {
            std::cout << "Need an int value for number of gens" << std::endl;
            return 1;
        }
        
        Shrinker *shrinker = new Shrinker(labelMapPath, labelNum, jacobianValue,
                padding, outDir, numGens);

}



