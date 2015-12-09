#include "itkImage.h"
#include "itkRandomImageSource.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
 
typedef itk::Image<float, 3> ImageType;


float getTextureFeature(ImageType::RegionType roi, ImageType::Pointer image)
{
    // 3 is inertia
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<ImageType> TextureFilterType;
  TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  textureFilter->SetNumberOfBinsPerAxis(16);
  textureFilter->SetPixelValueMinMax(0, 255);
  textureFilter->FastCalculationsOn();
  filter->SetRegionOfInterest(roi);
  filter->SetInput(image);
  textureFilter->SetInput(filter->GetOutput());
  textureFilter->Update();
 
  const TextureFilterType::FeatureValueVector* output = textureFilter->GetFeatureMeans();
  // defaults to: {Energy, Entropy, InverseDifferenceMoment, Inertia, ClusterShade, ClusterProminence
  
   
  return (*output)[3];
}


 
int main(int argc, char *argv[])
{
    if(argc < 2)
    {
    std::cerr << "Usage: " << argv[0] << " Required image.nrrd" << std::endl;
    return EXIT_FAILURE;
    }
    
    int chunkSize = 5;
    
    std::string fileName = argv[1];
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader=ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();
    ImageType::Pointer image=reader->GetOutput();
    //std::cout << (*output)[i] << std::endl;
            
    // Loop over image in chunks 
        //slide window over the entire image
   //image->GetLargestPossibleRegion().GetSize(0)-1; x++)
    
   
   ImageType::IndexType start;
   ImageType::SizeType size;
   size[0] = 5;
   size[1] = 5;
   size[2] = 5;
   
   float textureFeature;
              
   ImageType::SizeType imsize = image->GetLargestPossibleRegion().GetSize();
   std::cout << imsize << std::endl;
   
   for(int x=0; x<imsize[0] - chunkSize; x += chunkSize) // for all Columns
    {
        for(int y=0; y<imsize[1] - chunkSize; y+=chunkSize) // for all Rows
        {
            
            for(int z=0; z<imsize[2]- chunkSize; z+=chunkSize) // for all Rows
            {
               start[0] = x;
               start[1] = y;
               start[2] = z;
               //std::cout << x << ':' << y << ':' << z << std::endl;
               ImageType::RegionType desiredRegion;
               desiredRegion.SetSize(size);
               desiredRegion.SetIndex(start);
               
               //std::cout << 'index: ' << x << y << z;
               
               textureFeature = getTextureFeature(desiredRegion, image); 
               std::cout << textureFeature;
            }
        }
        
      }
    
        
}


