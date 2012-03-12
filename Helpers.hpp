#include "itkImageRegionConstIteratorWithIndex.h"

// STL
#include <algorithm>
#include <vector>
#include <numeric> // for accumulate()

// ITK
#include "itkImageFileWriter.h"

namespace Helpers
{
  
/** Copy the input to the output*/
template<typename TImage>
void DeepCopy(const TImage* const input, TImage* const output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

template<typename TPixel>
void DeepCopy(const itk::VectorImage<TPixel, 2>* const input,
              itk::VectorImage<TPixel, 2>* const output)
{
  typedef itk::VectorImage<TPixel, 2> ImageType;
  
  output->SetRegions(input->GetLargestPossibleRegion());
  output->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
  output->Allocate();

  itk::ImageRegionConstIterator<ImageType> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

template<typename TImage>
void WriteImage(const TImage* const image, const std::string& fileName)
{
  typedef  itk::ImageFileWriter< TImage > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName);
  writer->SetInput(image);
  writer->Update();
}

template<typename T>
void WriteVectorToFile(const std::vector<T> &v, const std::string& filename)
{
  std::ofstream fout(filename.c_str());
 
  for(unsigned int i = 0; i < v.size(); ++i)
    {
    fout << v[i] << std::endl;
    }
 
  fout.close();
}

template<typename T>
T VectorMedian(const std::vector<T> &v)
{
  
  int n = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
  

  
//   std::sort(v.begin(), v.end());
//   int n = v.size() / 2;
//   return v[n];
}

template<typename T>
T VectorAverage(const std::vector<T> &v)
{
  T vecSum = std::accumulate(v.begin(), v.end(), 0);
  
  return vecSum / static_cast<T>(v.size());
}

#if 0
// Specializations (must go before the calls to these functions - prevents "specialization after instantiation" errors)
// Specializationso RGBDI images are displayed using RGB only
template <>
void ITKImagetoVTKImage<RGBDIImageType>(RGBDIImageType::Pointer image, vtkImageData* outputImage)
{
  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(3); // we are definitly making an RGB image
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<RGBDIImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < 3; component++) // we explicitly stop at 3
      {
      pixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
      }

    ++imageIterator;
    }
}
#endif

template <typename TImage>
void SetPixels(TImage* const image, const std::vector<itk::Index<2> >& pixels, typename TImage::PixelType& value)
{
  for(unsigned int i = 0; i < pixels.size(); ++i)
    {
    image->SetPixel(pixels[i], value);
    }  
}

template <typename TImage>
void SetPixelsInRegionToValue(TImage* const image, const itk::ImageRegion<2>& region,
                              const typename TImage::PixelType& value)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(value);
    ++imageIterator;
    }
}


} // end namespace