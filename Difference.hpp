#ifndef Difference_HPP
#define Difference_HPP

#include <cmath>

class Difference
{
public:
  typedef itk::VariableLengthVector<float> VectorType;
  virtual float ComputeDifference(const VectorType& a,  const VectorType& b) = 0;
};

class DepthDifference : public Difference
{
  public:
  float ComputeDifference(const VectorType& a, const VectorType& b)
  {
    //std::cout << "Difference between " << a << " and " << b << std::endl;
    return pow(a[3] - b[3], 2);
  }
};

class ColorDifference : public Difference
{
  public:
  float ComputeDifference(const VectorType& a, const VectorType& b)
  {
    float sum = 0.0f;
    for(unsigned int component = 0; component < 3; ++component)
      {
      sum += pow(a[component] - b[component], 2);
      }
    return sum;
  }
};

class WeightedDifference : public Difference
{
  public:
  std::vector<float> Weights;
  
  WeightedDifference(const std::vector<float>& weights) : Weights(weights)
  {
  }
  
  float ComputeDifference(const VectorType& a, const VectorType& b)
  {
    float sum = 0.0f;
    for(unsigned int component = 0; component < a.GetSize(); ++component)
      {
      sum += Weights[component] * pow(a[component] - b[component], 2);
      }
    return sum;
  }
};

#endif
