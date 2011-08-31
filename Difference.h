#ifndef DIFFERENCE_H
#define DIFFERENCE_H

// STL
#include <cmath>

// VTK
#include <vtkMath.h>

// Custom
#include "Types.h"

class Difference
{
public:
  // Default constructor
  Difference();
  
  // Copy constructor
  Difference(const Difference& copy_from_me)
  {
    this->AverageColorDifference = copy_from_me.AverageColorDifference;
    this->MedianColorDifference = copy_from_me.MedianColorDifference;
    this->MinColorDifference = copy_from_me.MinColorDifference;
    this->MaxColorDifference = copy_from_me.MaxColorDifference;
  
    this->AverageDepthDifference = copy_from_me.AverageDepthDifference;
    this->MedianDepthDifference = copy_from_me.MedianDepthDifference;
    this->MinDepthDifference = copy_from_me.MinDepthDifference;
    this->MaxDepthDifference = copy_from_me.MaxDepthDifference;
    
  };
  
  
  void SetImage(ImageType::Pointer);
  
  void ComputeMinAndMaxInAllChannels();
  PixelType MinimumOfChannels;
  PixelType MaximumOfChannels;
  
  void WriteImages();
  
  void Compute();
  virtual float GetDifference(itk::Index<2>) = 0;
  //virtual float GetAverageDifference() = 0;
  float GetAverageDifference();
  
  FloatScalarImageType::Pointer NormalizedDepthDifferenceImage;
  FloatScalarImageType::Pointer NormalizedColorDifferenceImage;
  
  void CreateNormalizedDepthGradientMagnitude();
  void CreateNormalizedColorGradientMagnitude();
    
  float AverageColorDifference;
  float MedianColorDifference;
  float MinColorDifference;
  float MaxColorDifference;
  
  float AverageDepthDifference;
  float MedianDepthDifference;
  float MinDepthDifference;
  float MaxDepthDifference;
  
  float AverageDifference;
  
protected:
  void CreateRGBImage(Vector3ImageType::Pointer image);
  void CreateDepthImage(FloatScalarImageType::Pointer image);
  
  ImageType::Pointer Image;
};

class DifferenceDepth : public Difference
{
public:
  // Copy constructor
  DifferenceDepth(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceDepth(){}
  
  float GetDifference(itk::Index<2> pixel)
  {
    return this->NormalizedDepthDifferenceImage->GetPixel(pixel);
  }
};


class DifferenceColor : public Difference
{
public:
  // Copy constructor
  DifferenceColor(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceColor(){}
  
  float GetDifference(itk::Index<2> pixel)
  {
    return this->NormalizedColorDifferenceImage->GetPixel(pixel);
  }
};


class DifferenceMaxOfColorOrDepth : public Difference
{
public:
  // Copy constructor
  DifferenceMaxOfColorOrDepth(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceMaxOfColorOrDepth(){}
  
  float GetDifference(itk::Index<2> pixel)
  {
    return std::max(this->NormalizedColorDifferenceImage->GetPixel(pixel), this->NormalizedDepthDifferenceImage->GetPixel(pixel));
  }
  
  
};
  
#if 0
class DifferenceDepthWeightedByColor : public Difference
{
public:
  // Copy constructor
  DifferenceDepthWeightedByColor(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceDepthWeightedByColor(){}
  
  float Compute(PixelType a, PixelType b)
  {
    // Using normalized values is much better because then the main graph cut lambda does not change wildly from data set to data set
    DifferenceColorDataNormalized differenceColorDataNormalized(*this);
    float colorDifference = differenceColorDataNormalized.Compute(a,b);
    
    DifferenceDepthDataNormalized differenceDepthDataNormalized(*this);
    float depthDifference = differenceDepthDataNormalized.Compute(a,b);

    return std::max(depthDifference, colorDifference) + (depthDifference + colorDifference)/2.0;
  }
};

class DifferenceEuclidean : public Difference
{
public:
  // Copy constructor
  DifferenceEuclidean(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceEuclidean(){}

  float Compute(PixelType a, PixelType b)
  {
    float difference = 0.0;
    for(unsigned int i = 0; i < 3u; i++)
      {
      difference += pow(a[i] - b[i],2);
      }
    return sqrt(difference);
  }

};
#endif

#endif
