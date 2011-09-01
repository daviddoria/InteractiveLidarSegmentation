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
unsigned int CountNonZeroPixels(MaskImageType::Pointer image);
  
bool FindClosestNonZeroPixel(MaskImageType::Pointer image, itk::Index<2> queryPixel, unsigned int radiusValue, itk::Index<2>&);
itk::Index<2> FindClosestNonZeroPixel(MaskImageType::Pointer, itk::Index<2>);

void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage);

void CreateTransparentImage(vtkImageData* VTKImage);

void SetImageSize(vtkImageData* input, vtkImageData* output);

void SetPixels(vtkImageData* VTKImage, std::vector<itk::Index<2> > pixels, unsigned char color[3]);

// Determine if a number is NaN
bool IsNaN(const double a);

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image);

// Create a list of the non-zero pixels in 'image'
std::vector<itk::Index<2> > BinaryImageToIndices(UnsignedCharScalarImageType::Pointer image);

// Invert binary image
void InvertBinaryImage(UnsignedCharScalarImageType::Pointer image, UnsignedCharScalarImageType::Pointer inverted);

void ITKImagetoVTKImage(ImageType::Pointer image, vtkImageData* outputImage); // This function simply drives ITKImagetoVTKRGBImage or ITKImagetoVTKMagnitudeImage
void ITKImagetoVTKRGBImage(ImageType::Pointer image, vtkImageData* outputImage);
void ITKImagetoVTKMagnitudeImage(ImageType::Pointer image, vtkImageData* outputImage);

void ITKScalarImageToVTKImage(MaskImageType::Pointer image, vtkImageData* outputImage);


std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata);

template<typename T>
T VectorMedian(std::vector<T> &v);

template<typename T>
T VectorAverage(std::vector<T> &v);

template<typename TImage>
void WriteImage(typename TImage::Pointer image, const std::string& fileName);

}

#include "Helpers.txx"

#endif
