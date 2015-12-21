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

typedef itk::Image<float, 3> FloatImageType;

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<VectorType, 3> VectorImageType;


typedef itk::DisplacementFieldJacobianDeterminantFilter<VectorImageType >  JacFilterType;
typedef itk::SquaredDifferenceImageFilter<FloatImageType, FloatImageType, FloatImageType> SquareDiffFilter;
typedef itk::StatisticsImageFilter<FloatImageType> StatsFilter;

/*
 * 
 */

const float MUTATE_VAL = 30.0;

class Shrinker
{
    
public:
    FloatImageType::Pointer idealJacImage; 
    VectorImageType::Pointer defField;
    Shrinker(string idealPath, string outDir, int numGens);   
    void run();
    float calculateFitness();
    void mutate();
    void makeDefField();
    float randomFloat(float, float);
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

void Shrinker::run(){
    
    this->makeDefField();

    int i = 0;
    while(i < this->numGens){
       this->mutate();
       i++;
       if (i % 100 == 0){
           std::cout << "done: " << std::cout << i << std::endl;
       }
    };
    
}

void Shrinker::makeDefField(){
    // Make a zero vector field same size as the ideal jacobian image
    FloatImageType::SizeType jacSize;
    jacSize = this->idealJacImage->GetLargestPossibleRegion().GetSize();
    VectorImageType::Pointer def = VectorImageType::New();
    VectorImageType::SizeType defSize;
    
    VectorImageType::RegionType region;
    
    defSize[0] = jacSize[0];
    defSize[1] = jacSize[1];
    defSize[2] = jacSize[2];
    
    
    region.SetSize(defSize);
    
    def->SetRegions(region);
    def->Allocate();
    this->defField = def;
}

float Shrinker::calculateFitness(){
    
    JacFilterType::Pointer jacFilter = JacFilterType::New();
    jacFilter->SetInput(this->defField);
    FloatImageType::Pointer generatedJac = FloatImageType::New();
    //generatedJac = jacFilter->GetOutput();
    
    SquareDiffFilter::Pointer sqFilter = SquareDiffFilter::New();
    sqFilter->SetInput1(jacFilter->GetOutput());
    sqFilter->SetInput2(this->idealJacImage);
    
    FloatImageType::Pointer sqOut;
    sqOut = sqFilter->GetOutput();
    
    StatsFilter::Pointer statsFilter = StatsFilter::New();
    statsFilter->SetInput(sqOut);
    float sum = statsFilter->GetSum();
  
    return sum;
}

void Shrinker::mutate(){
    // Mutate a random vector element
    while (true){
        int xRandPos = rand() % this->xShape;
        int yRandPos = rand() % this->yShape;
        int zRandPos = rand() % this->zShape;
        int vectorRandPos = rand() % 3;

        VectorImageType::IndexType vidx;
        vidx[0] = xRandPos;
        vidx[1] = yRandPos;
        vidx[2] = zRandPos;

        VectorImageType::PixelType original_vec;

        float fitness = this->fitness;
    
        // Get the original def component value
        original_vec = this->defField->GetPixel(vidx);
        float original_val = original_vec[vectorRandPos];
    
       
       float mutateRange = this->randomFloat(0.0, this->fitness * 30.0);
       float randVal = this->randomFloat(-mutateRange, mutateRange);
       original_vec[vectorRandPos] = randVal;
       this->defField->SetPixel(vidx, original_vec);
       
       float newFitness = this->calculateFitness();
       std::cout << this->fitness << std::endl;
       if (newFitness < this->fitness){
           this->fitness == newFitness;
           break;
       }
       else{
           // replace the original value
           original_vec[vectorRandPos] = original_val;
           this->defField->SetPixel(vidx, original_val);
       }
       
    }
    this->fitness = fitness;
    return; 
}


float Shrinker::randomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
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



