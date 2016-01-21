/* 
 * File:   main.cpp
 * Author: neil
 *
 * Created on 17 December 2015, 12:00
 */
#include "itkImage.h"
#include "itkRandomImageSource.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream  
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>

using namespace std;


typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<VectorType, 3> VectorImageType;
typedef itk::DisplacementFieldJacobianDeterminantFilter<VectorImageType >  JacFilterType;
typedef itk::Image<float, 3> FloatImageType;



int main(int argc, char** argv) {
    
    std::string idealJacPath = argv[1];
    std::string numGensString = argv[2];
    
    int numGens;
    istringstream ( numGensString ) >> numGens;

	
        
    typedef itk::ImageFileReader<VectorImageType> VectorReaderType;
    VectorReaderType::Pointer reader=VectorReaderType::New();
    reader->SetFileName(idealJacPath);
    reader->Update();
    VectorImageType::Pointer vectorImage=reader->GetOutput();
    
    
    JacFilterType::Pointer jacFilter = JacFilterType::New();

    FloatImageType::Pointer generatedJac = FloatImageType::New();
    
    for (int i =0; i < numGens; i++){
        jacFilter->SetInput(vectorImage);
        jacFilter->Update();
        jacFilter->Modified();
        generatedJac = jacFilter->GetOutput();
    }
    return 0;

}



