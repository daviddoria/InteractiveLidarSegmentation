#ifndef DIFFERENCEMANUAL_H
#define DIFFERENCEMANUAL_H

#include <cmath>

#include <vtkMath.h>

class CIELABColorDifference;
class HSVColorDifference;
class DifferenceDepth;
class DifferenceColor;
class DifferenceMaxOfColorOrDepth;
class DifferenceDepthWeightedByColor;
class DifferenceColorDataNormalized;
class DifferenceColorAbsoluteNormalized;
class DifferenceDepthDataNormalized;
class DifferenceDepthAbsoluteNormalized;

class Difference
{
public:
  // Default constructor
  Difference(){};
  
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
  
  virtual float Compute(PixelType a, PixelType b) = 0;
  
  float AverageColorDifference;
  float MedianColorDifference;
  float MinColorDifference;
  float MaxColorDifference;
  
  float AverageDepthDifference;
  float MedianDepthDifference;
  float MinDepthDifference;
  float MaxDepthDifference;
  
protected:

};

class DifferenceDepth : public Difference
{
public:
  // Copy constructor
  DifferenceDepth(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceDepth(){}
  
  float Compute(PixelType a, PixelType b)
  {
    // Compute the Euclidean distance between N dimensional pixels
    //float difference = pow(a[3] - b[3],2);
    //return sqrt(difference);

    float difference = fabs(a[3] - b[3]);
    //difference = std::min(difference, 10.0f); // Invalid returns sometimes cause erroneously high depth differences, so clamp them to [0,10]
      
    return difference;
  }
};


class DifferenceColor : public Difference
{
public:
  // Copy constructor
  DifferenceColor(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceColor(){}
  
  float Compute(PixelType a, PixelType b)
  {
    float difference = 0;

    // Compute the Euclidean RGB distance between N dimensional pixels
    //for(unsigned int i = 0; i < std::min(this->Image->GetNumberOfComponentsPerPixel(), 3u); i++)
    for(unsigned int i = 0; i < 3u; i++)
      {
      difference += fabs(a[i] - b[i]);
      }

    return difference;
  }
};


class DifferenceColorDataNormalized : public Difference
{
public:
  // Copy constructor
  DifferenceColorDataNormalized(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceColorDataNormalized(){}
  
  float Compute(PixelType a, PixelType b)
  {
    DifferenceColor differenceFunction;
    float difference = differenceFunction.Compute(a,b);
    float normalizedDifference = (difference - this->MinColorDifference)/(this->MaxColorDifference - this->MinColorDifference);
    return normalizedDifference;
  }
};

class DifferenceColorAbsoluteNormalized : public Difference
{
public:
  // Copy constructor
  DifferenceColorAbsoluteNormalized(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceColorAbsoluteNormalized(){}
  
  float Compute(PixelType a, PixelType b)
  {
    DifferenceColor differenceFunction;
    float difference = differenceFunction.Compute(a,b);
    float normalizedDifference = difference/(255.*3.);
    return normalizedDifference;
  }
};

class DifferenceDepthDataNormalized : public Difference
{
public:
  // Copy constructor
  DifferenceDepthDataNormalized(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceDepthDataNormalized(){}
  
  float Compute(PixelType a, PixelType b)
  {
    DifferenceDepth differenceFunction;
    float difference = differenceFunction.Compute(a,b);
    
    float normalizedDifference = (difference - this->MinDepthDifference)/(this->MaxDepthDifference - this->MinDepthDifference);

    return normalizedDifference;
  }
};

class DifferenceDepthAbsoluteNormalized : public Difference
{
public:
  // Copy constructor
  DifferenceDepthAbsoluteNormalized(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceDepthAbsoluteNormalized(){}
  
  float Compute(PixelType a, PixelType b)
  {
    DifferenceDepth differenceFunction;
    float difference = differenceFunction.Compute(a,b);
    
    float normalizedDifference = difference/10.0; // Normalize to the arbitrarily decided maximum possible depth difference
    
    return normalizedDifference;
  }
};

class HSVColorDifference : public Difference
{
public:
  // Copy constructor
  HSVColorDifference(const Difference& input) : Difference(input) {}
  
  // Default constructor
  HSVColorDifference(){}
  
  float Compute(PixelType a, PixelType b)
  {
    float rgbA[3];
    rgbA[0] = a[0];
    rgbA[1] = a[1];
    rgbA[2] = a[2];

    float rgbB[3];
    rgbB[0] = b[0];
    rgbB[1] = b[1];
    rgbB[2] = b[2];

    // Convert to HSV
    float hsvA[3];
    vtkMath::RGBToHSV(rgbA, hsvA);
    
    float hsvB[3];
    vtkMath::RGBToHSV(rgbB, hsvB);
    
    float difference = 0.0;
    // Compute the HSV Euclidean distance
    for(unsigned int i = 0; i < 3u; i++)
      {
      difference += fabs(hsvA[i] - hsvB[i]);
      }
    
    // Compute the hue Euclidean distance
    //difference += fabs(hsvA[0] - hsvB[0]);
    
    return difference;
  }
};

class CIELABColorDifference : public Difference
{
public:
  // Copy constructor
  CIELABColorDifference(const Difference& input) : Difference(input) {}
  
  // Default constructor
  CIELABColorDifference(){}
  
  float Compute(PixelType a, PixelType b)
  {
    double rgbA[3];
    rgbA[0] = a[0];
    rgbA[1] = a[1];
    rgbA[2] = a[2];

    double rgbB[3];
    rgbB[0] = b[0];
    rgbB[1] = b[1];
    rgbB[2] = b[2];

    // Convert to CIELAB
    double labA[3];
    vtkMath::RGBToLab(rgbA, labA);
    
    double labB[3];
    vtkMath::RGBToLab(rgbB, labB);
    
    float difference = 0.0;
    //Compute CIELAB Euclidean distance
    for(unsigned int i = 0; i < 3u; i++)
      {
      difference += fabs(labA[i] - labB[i]);
      }
    
    // The a* and b* (elements 1 and 2) encode color. Element 0 (not used here), L*, encodes lightness.
    //difference = fabs(labA[1] - labB[1]) + fabs(labA[2] - labB[2]);
    
    return difference;
  }
};


class DifferenceMaxOfColorOrDepth : public Difference
{
public:
  // Copy constructor
  DifferenceMaxOfColorOrDepth(const Difference& input) : Difference(input) {}
  
  // Default constructor
  DifferenceMaxOfColorOrDepth(){}
  
  float Compute(PixelType a, PixelType b)
  {
    DifferenceColorDataNormalized differenceColorDataNormalized;
    float colorDifference = differenceColorDataNormalized.Compute(a,b);
    
    DifferenceDepthDataNormalized differenceDepthDataNormalized;
    float depthDifference = differenceDepthDataNormalized.Compute(a,b);

    return std::max(colorDifference, depthDifference);
  }
};
  
  
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
