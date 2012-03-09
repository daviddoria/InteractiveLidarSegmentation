#ifndef HELPERS_H
#define HELPERS_H

// ITK
#include "itkImage.h"
#include "itkIndex.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

// Custom
#include "Types.h"

class vtkPolyData;

namespace Helpers
{
unsigned int CountNonZeroPixels(const MaskImageType* const image);

float NegativeLog(const float);

std::vector<itk::Index<2> > GetNonZeroPixels(const MaskImageType* const image);
  
bool FindClosestNonZeroPixel(const MaskImageType* const image, itk::Index<2> queryPixel, const unsigned int radiusValue, const itk::Index<2>& index);
itk::Index<2> FindClosestNonZeroPixel(const MaskImageType* const mask, itk::Index<2>& returnIndex);

void MaskImage(vtkImageData* const VTKImage, vtkImageData* const VTKSegmentMask, vtkImageData* const VTKMaskedImage);

void CreateTransparentImage(vtkImageData* const VTKImage);

void SetImageSize(vtkImageData* const input, vtkImageData* const output);

void SetPixels(vtkImageData* const VTKImage, const std::vector<itk::Index<2> >& pixels, const unsigned char color[3]);

// Determine if a number is NaN
bool IsNaN(const double a);

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(const std::vector<itk::Index<2> >& indices, UnsignedCharScalarImageType* const image);

// Create a list of the non-zero pixels in 'image'
std::vector<itk::Index<2> > BinaryImageToIndices(const UnsignedCharScalarImageType* const image);

// Invert binary image
void InvertBinaryImage(const UnsignedCharScalarImageType* const image, UnsignedCharScalarImageType* const inverted);

void ITKImagetoVTKImage(const ImageType* const image, vtkImageData* const outputImage); // This function simply drives ITKImagetoVTKRGBImage or ITKImagetoVTKMagnitudeImage
void ITKImagetoVTKRGBImage(const ImageType* const image, vtkImageData* const outputImage);
void ITKImagetoVTKMagnitudeImage(const ImageType* const image, vtkImageData* const outputImage);

void ITKScalarImageToVTKImage(const MaskImageType* const image, vtkImageData* const outputImage);


std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* const polydata);

template<typename T>
T VectorMedian(const std::vector<T> &v);

template<typename T>
T VectorAverage(const std::vector<T> &v);

template<typename TImage>
void WriteImage(const TImage* const image, const std::string& fileName);

template <typename TImage>
void SetPixels(TImage* image, const std::vector<itk::Index<2> >& pixels, typename TImage::PixelType& value);

template<typename TImage>
void DeepCopy(const TImage* const input, TImage* const output);

template<typename TImage>
void DeepCopyVectorImage(const TImage* const input, TImage* const output);
}

#include "Helpers.txx"

#endif
