#include "Difference.h"
#include "Helpers.h"

// ITK
#include "itkGradientMagnitudeImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

Difference::Difference()
{
  this->NormalizedDepthDifferenceImage = NULL;
  this->NormalizedColorDifferenceImage = NULL;
  
  //this->Image = ImageType::New();
  
}

float Difference::GetAverageDifference()
{
  itk::ImageRegionIterator<ImageType> imageIterator(this->Image, this->Image->GetLargestPossibleRegion());
 
  float sumOfDifferences = 0.;
  unsigned int numberOfDifferences = 0;
  
  while(!imageIterator.IsAtEnd())
    {
    sumOfDifferences += this->GetDifference(imageIterator.GetIndex());
    numberOfDifferences++;
    ++imageIterator;
    }
  return sumOfDifferences/static_cast<float>(numberOfDifferences);
}

void Difference::SetImage(ImageType::Pointer image)
{
  this->Image = image;
  //Helpers::DeepCopy<ImageType>(image, this->Image);
  
  NormalizedDepthDifferenceImage = FloatScalarImageType::New();  
  NormalizedColorDifferenceImage = FloatScalarImageType::New();
  
  Compute();
  
  this->AverageDifference = GetAverageDifference();
}

void Difference::Compute()
{
  if(!this->Image)
    {
    throw std::runtime_error("An image must be set before calling Difference::Compute()!");
    }

  CreateNormalizedColorGradientMagnitude();

  CreateNormalizedDepthGradientMagnitude();
  ComputeMinAndMaxInAllChannels();  
}

void Difference::ComputeMinAndMaxInAllChannels()
{
  std::cout << "ComputeMinAndMaxInAllChannels()" << std::endl;
  
  this->MinimumOfChannels.SetSize(this->Image->GetNumberOfComponentsPerPixel());
  this->MaximumOfChannels.SetSize(this->Image->GetNumberOfComponentsPerPixel());
  
  for(unsigned int i = 0; i < this->Image->GetNumberOfComponentsPerPixel(); ++i)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<ImageType, FloatScalarImageType> IndexSelectionType;
    IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(this->Image);
    indexSelectionFilter->Update();
  
    typedef itk::MinimumMaximumImageCalculator <FloatScalarImageType> ImageCalculatorFilterType;
    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
    imageCalculatorFilter->SetImage(indexSelectionFilter->GetOutput());
    imageCalculatorFilter->Compute();
    
    this->MinimumOfChannels[i] = imageCalculatorFilter->GetMinimum();
    this->MaximumOfChannels[i] = imageCalculatorFilter->GetMaximum();
    }
}


void Difference::WriteImages()
{
  Helpers::WriteImage(NormalizedColorDifferenceImage.GetPointer(), "colorGradientMagnitude.mha");
  
  // Compute depth gradient
  Helpers::WriteImage(NormalizedDepthDifferenceImage.GetPointer(), "depthGradientMagnitude.mha");
  
  typedef itk::MaximumImageFilter<FloatScalarImageType> MaximumImageFilterType;
  MaximumImageFilterType::Pointer maximumImageFilter = MaximumImageFilterType::New ();
  maximumImageFilter->SetInput(0, NormalizedColorDifferenceImage);
  maximumImageFilter->SetInput(1, NormalizedDepthDifferenceImage);
  maximumImageFilter->Update();
  Helpers::WriteImage<FloatScalarImageType>(maximumImageFilter->GetOutput(), "maxGradientMagnitude.mha");
  
#if 0
  // Combine cleverly
  FloatScalarImageType::Pointer combinedImage = FloatScalarImageType::New();
  combinedImage->SetRegions(this->Image->GetLargestPossibleRegion());
  combinedImage->Allocate();
  combinedImage->FillBuffer(0);
  
  itk::ImageRegionIterator<FloatScalarImageType> combinedImageIterator(combinedImage, combinedImage->GetLargestPossibleRegion());
 
  while(!combinedImageIterator.IsAtEnd())
    {
    
//     float depthDifference = depthGradientMagnitudeImage->GetPixel(combinedImageIterator.GetIndex());
//     float colorDifference = rgbGradientMagnitudeImage->GetPixel(combinedImageIterator.GetIndex());
//     
//     //float pixelDifference = std::max(depthDifference, colorDifference) + (depthDifference + colorDifference)/2.0;
//     float pixelDifference = 0.0;
//     if(depthDifference < .02)
//       {
//       pixelDifference = depthDifference; // If the depth difference is very small then we are in a flat region (??? this is not true, we could also be at a ground boundary)
//       }
//     else
//       {
//       pixelDifference = depthDifference;
//       }
//       
//     combinedImageIterator.Set(pixelDifference);
//  
    ++combinedImageIterator;
    }
  Helpers::WriteImage<FloatScalarImageType>(combinedImage, "combinedGradientMagnitude.mha");
  
#endif
}


void Difference::CreateRGBImage(Vector3ImageType::Pointer image)
{
  image->SetRegions(this->Image->GetLargestPossibleRegion());
  image->Allocate();
  
  Vector3ImageType::PixelType zeroPixel;
  zeroPixel.Fill(0);
  
  image->FillBuffer(zeroPixel);
  
  itk::ImageRegionIterator<ImageType> fullImageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Vector3ImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  
  while(!fullImageIterator.IsAtEnd())
    {
    ImageType::PixelType fullPixel = fullImageIterator.Get();
  
    Vector3ImageType::PixelType rgbPixel;
    for(unsigned int i = 0; i < 3; ++i)
      {
      rgbPixel[i] = fullPixel[i];
      }
  
    rgbImageIterator.Set(rgbPixel);
    ++fullImageIterator;
    ++rgbImageIterator;
    }
}


void Difference::CreateDepthImage(FloatScalarImageType::Pointer image)
{
  image->SetRegions(this->Image->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);
  
  itk::ImageRegionIterator<ImageType> fullImageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatScalarImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  
  while(!fullImageIterator.IsAtEnd())
    {
    ImageType::PixelType fullPixel = fullImageIterator.Get();
  
    float depthPixel = fullPixel[3];
  
    rgbImageIterator.Set(depthPixel);
    ++fullImageIterator;
    ++rgbImageIterator;
    }
}

void Difference::CreateNormalizedColorGradientMagnitude()
{
  Vector3ImageType::Pointer rgbImage = Vector3ImageType::New();
  CreateRGBImage(rgbImage);
  
//   typedef itk::BilateralImageFilter<Vector3ImageType, Vector3ImageType> BilateralFilterType;
//   BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
//   bilateralFilter->SetInput( rgbImage);
//   bilateralFilter->SetDomainSigma(4);
//   bilateralFilter->SetRangeSigma(10);
//   bilateralFilter->Update();
  
  typedef itk::Image< itk::CovariantVector<float, 3>,  2 >    VectorImageType;
  typedef itk::VectorGradientMagnitudeImageFilter<VectorImageType>  ColorVectorGradientMagnitudeImageFilterType;
  ColorVectorGradientMagnitudeImageFilterType::Pointer colorGradientFilter = ColorVectorGradientMagnitudeImageFilterType::New();
  colorGradientFilter->SetInput(rgbImage);
  //colorGradientFilter->SetInput(bilateralFilter->GetOutput());
  colorGradientFilter->SetUsePrincipleComponentsOff();
  colorGradientFilter->Update();
  
  typedef itk::RescaleIntensityImageFilter< ColorVectorGradientMagnitudeImageFilterType::OutputImageType, FloatScalarImageType > ColorGradientRescaleFilterType;
  ColorGradientRescaleFilterType::Pointer colorGradientRescaleFilter = ColorGradientRescaleFilterType::New();
  colorGradientRescaleFilter->SetInput(colorGradientFilter->GetOutput());
  colorGradientRescaleFilter->SetOutputMinimum(0);
  colorGradientRescaleFilter->SetOutputMaximum(1);
  colorGradientRescaleFilter->Update();
  
  Helpers::DeepCopy<FloatScalarImageType>(colorGradientRescaleFilter->GetOutput(), NormalizedColorDifferenceImage);
}

void Difference::CreateNormalizedDepthGradientMagnitude()
{
//   FloatScalarImageType::Pointer depthImage = FloatScalarImageType::New();
//   CreateDepthImage(depthImage);
//   
//   typedef itk::GradientMagnitudeImageFilter<FloatScalarImageType, FloatScalarImageType >  DepthGradientMagnitudeImageFilterType;
//   DepthGradientMagnitudeImageFilterType::Pointer depthGradientMagnitudeImageFilter = DepthGradientMagnitudeImageFilterType::New();
//   depthGradientMagnitudeImageFilter->SetInput(depthImage);
//   depthGradientMagnitudeImageFilter->Update();
//   
//   typedef itk::RescaleIntensityImageFilter< FloatScalarImageType, FloatScalarImageType > DepthGradientRescaleFilterType;
//   DepthGradientRescaleFilterType::Pointer depthGradientRescaleFilter = DepthGradientRescaleFilterType::New();
//   depthGradientRescaleFilter->SetInput(depthGradientMagnitudeImageFilter->GetOutput());
//   depthGradientRescaleFilter->SetOutputMinimum(0);
//   depthGradientRescaleFilter->SetOutputMaximum(1);
//   depthGradientRescaleFilter->Update();
//   
//   Helpers::DeepCopy(depthGradientRescaleFilter->GetOutput(), NormalizedDepthDifferenceImage.GetPointer());
}
