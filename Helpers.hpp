#include "itkImageRegionConstIteratorWithIndex.h"

// STL
#include <algorithm>
#include <vector>
#include <numeric> // for accumulate()

// ITK
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

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

template<typename TImage>
void ITKScalarImageToVTKImage(const TImage* const image, vtkImageData* const outputImage)
{
  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<
                  TImage, UnsignedCharScalarImageType > RescaleFilterType;

  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(image);
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType>
    imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }
}


template<typename TPixel>
void ExtractChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                    itk::Image<TPixel, 2>* const output)
{
  typedef itk::VectorImage<TPixel, 2> VectorImageType;
  typedef itk::Image<TPixel, 2> ScalarImageType;

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(channel);
  indexSelectionFilter->SetInput(image);
  indexSelectionFilter->Update();

  DeepCopy(indexSelectionFilter->GetOutput(), output);
}

template<typename TImage>
void NormalizeImage(const TImage* const image, TImage* const outputImage)
{
  DeepCopy(image, outputImage);
  
  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    std::cout << "Normalizing channel " << channel << std::endl;
    typedef itk::Image<float, 2> ScalarImageType;
    ScalarImageType::Pointer scalarImage = ScalarImageType::New();
    ExtractChannel(image, channel, scalarImage.GetPointer());
    
    float mean = MeanValue(scalarImage.GetPointer());
    std::cout << "Channel " << channel << " mean is " << mean << std::endl;
    
//     float variance = Variance(scalarImage.GetPointer());
//     std::cout << "Channel " << channel << " variance is " << variance << std::endl;
    float standardDeviation = StandardDeviation(scalarImage.GetPointer());
    std::cout << "Channel " << channel << " standardDeviation is " << standardDeviation << std::endl;
    
    itk::ImageRegionIterator<ScalarImageType> imageIterator(scalarImage, scalarImage->GetLargestPossibleRegion());
    while(!imageIterator.IsAtEnd())
      {
      float newValue = (imageIterator.Get() - mean) / standardDeviation;
      imageIterator.Set(newValue);
      ++imageIterator;
      }
    ReplaceChannel(outputImage, channel, scalarImage.GetPointer(), outputImage);
  }
}

// template<typename TPixel>
// void ReplaceChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
//                     const itk::Image<TPixel, 2>* const replacement, itk::VectorImage<TPixel, 2>* const output)
// {
//   if(image->GetLargestPossibleRegion() != replacement->GetLargestPossibleRegion())
//     {
//     throw std::runtime_error("Image and replacement channel are not the same size!");
//     }
// 
//   DeepCopy(image, output);
// 
//   itk::ImageRegionConstIterator<itk::VectorImage<TPixel, 2> > iterator(image, image->GetLargestPossibleRegion());
// 
//   while(!iterator.IsAtEnd())
//     {
//     typename itk::VectorImage<TPixel, 2>::PixelType pixel = iterator.Get();
//     pixel[channel] = replacement->GetPixel(iterator.GetIndex());
//     output->SetPixel(iterator.GetIndex(), pixel);
//     ++iterator;
//     }
// }

template<typename TImage>
float Variance(const TImage* const image)
{
  float average = MeanValue(image);

  float channelVarianceSummation = 0.0f;

  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    channelVarianceSummation += pow(imageIterator.Get() - average, 2);
    ++imageIterator;
    }
  float variance = channelVarianceSummation / static_cast<float>(image->GetLargestPossibleRegion().GetNumberOfPixels() - 1); // This (N-1) term in the denominator is for the "unbiased" sample variance. This is what is used by Matlab, Wolfram alpha, etc.

  return variance;
}

template<typename TImage>
float StandardDeviation(const TImage* const image)
{
  return sqrt(Variance(image));
}

template<typename TImage>
float MeanValue(const TImage* const image)
{
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  float sum = 0.0f;
  while(!imageIterator.IsAtEnd())
    {
    sum += imageIterator.Get();

    ++imageIterator;
    }
  return sum / static_cast<float>(image->GetLargestPossibleRegion().GetNumberOfPixels());
}

} // end namespace
