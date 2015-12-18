/* 
 * File:   main.cpp
 * Author: neil
 *
 * Created on 17 December 2015, 12:00
 */
#include "itkImage.h"
#include "itkRandomImageSource.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
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
typedef itk::Image<float, 4> VectorImageType;
typedef itk::DisplacementFieldJacobianDeterminantFilter<
               VectorImageType >  JacFilterType;

/*
 * 
 */


class Shrinker
{
    
public:
    FloatImageType::Pointer idealJacImage; 
    VectorImageType::Pointer defField;
    Shrinker(string idealPath, string outDir, int numGens);   
    void run();
    float calculateFitness();
    void mutate(FloatImageType::Pointer);
    void makeDefField();
    VectorImageType::IndexType getRandomPosition();
    int xShape;
    int yShape;
    int zShape;
    
    int numGens;
    float fitness;
    
};


Shrinker::Shrinker(string idealJacPath, string outDir, int numGens){
    
    typedef itk::ImageFileReader<FloatImageType> ByteReaderType;
    ByteReaderType::Pointer reader=ByteReaderType::New();
    reader->SetFileName(idealJacPath);
    reader->Update();
    FloatImageType::Pointer idealJacImage=reader->GetOutput();
    this->idealJacImage = idealJacImage;
    this->numGens = numGens;
    
    FloatImageType::SizeType jacSize = idealJacImage->GetLargestPossibleRegion().GetSize();
    xShape = jacSize[0];
    yShape = jacSize[1];
    zShape = jacSize[2];
 }

void Shrinker::makeDefField(){
    FloatImageType::SizeType jacSize;
    jacSize = this->idealJacImage->GetLargestPossibleRegion().GetSize();
    VectorImageType::Pointer def = VectorImageType::New();
    VectorImageType::SizeType defSize;
    
    VectorImageType::RegionType region;
    
    defSize[0] = jacSize[0];
    defSize[1] = jacSize[1];
    defSize[2] = jacSize[2];
    defSize[3] = 3;
    
    region.SetSize(defSize);
    
    def->SetRegions(region);
    def->Allocate();
    this->defField = def;
    
}

void Shrinker::run(){
    
    this->makeDefField();
    int i = 0;
    while(i < this->numGens){
        std::cout << "hello" << std::endl;
    };
    
}

float Shrinker::calculateFitness(){
    
    JacFilterType::Pointer filter = JacFilterType::New();
    filter->SetInput(this->defField);
    FloatImageType::Pointer generatedJac = FloatImageType::New();
    generatedJac = filter->GetOutput();
    float fitness = 4.0;
    return fitness;
}

void Shrinker::mutate(FloatImageType::Pointer jacImage){
    // Mutate a random vector element
    int xRand = rand() % this->xShape;
    int yRand = rand() % this->yShape;
    int zRand = rand() % this->zShape;
    int vectorRand = rand() % 3;
    
    float fitness = this->fitness;
    while (fitness >= this->fitness){
       fitness =  this->calculateFitness();
    }
    this->fitness = fitness;
    return;
    
}



int main(int argc, char** argv) {

    if(argc < 3)
    {
    std::cerr << "Usage:\n" << "ideal jac image path<string>\nnout folder<string>\nnumber of generations<int>\n" << std::endl;
    return EXIT_FAILURE;
    }
    
    std::string idealJacPath = argv[1];
    
    std::string outDir = argv[2];
     
    std::istringstream iss( argv[3] );
        int numGens;

        if (iss >> numGens)
        {
            std::cout << numGens << std::endl;
        }else
        {
            std::cout << "Need an int value for number of gens" << std::endl;
            return 1;
        }
        
        Shrinker *shrinker = new Shrinker(idealJacPath, outDir, numGens);
        shrinker->run();

}



