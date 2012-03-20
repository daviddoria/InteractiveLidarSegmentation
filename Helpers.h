#ifndef HELPERS_H
#define HELPERS_H

// ITK
#include "itkImage.h"
#include "itkIndex.h"
#include "itkVectorImage.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

// Custom
#include "Types.h"

class vtkPolyData;

namespace Helpers
{

/** We need a 1x1 radius quite often, so we can create it in a single line with this function */
itk::Size<2> Get1x1Radius();

std::vector<float> ComputeMinOfAllChannels(const ImageType* const image);
std::vector<float> ComputeMaxOfAllChannels(const ImageType* const image);

unsigned int CountNonZeroPixels(const MaskImageType* const image);

float NegativeLog(const float);

std::vector<itk::Index<2> > GetNonZeroPixels(const MaskImageType* const image);
  
bool FindClosestNonZeroPixel(const MaskImageType* const image, itk::Index<2> queryPixel,
                             const unsigned int radiusValue, const itk::Index<2>& index);
itk::Index<2> FindClosestNonZeroPixel(const MaskImageType* const mask, itk::Index<2>& returnIndex);

void MaskImage(vtkImageData* const VTKImage, vtkImageData* const VTKSegmentMask, vtkImageData* const VTKMaskedImage);

void MaskImage(const itk::VectorImage<float, 2>* const image, const MaskImageType* const mask, itk::VectorImage<float, 2>* const maskedImage);

void CreateTransparentImage(vtkImageData* const VTKImage);

void SetImageSize(vtkImageData* const input, vtkImageData* const output);

void SetPixels(vtkImageData* const VTKImage, const std::vector<itk::Index<2> >& pixels,
               const unsigned char color[3]);

/** Determine if a number is NaN */
bool IsNaN(const double a);

/** Mark each pixel at the specified 'indices' as a non-zero pixel in 'image' */
void IndicesToBinaryImage(const std::vector<itk::Index<2> >& indices, UnsignedCharScalarImageType* const image);

/** Create a list of the non-zero pixels in 'image' */
std::vector<itk::Index<2> > BinaryImageToIndices(const UnsignedCharScalarImageType* const image);

/** Invert a binary image */
void InvertBinaryImage(const UnsignedCharScalarImageType* const image, UnsignedCharScalarImageType* const inverted);

/** This function simply drives ITKImageToVTKRGBImage or ITKImageToVTKMagnitudeImage */
void ITKImageToVTKImage(const ImageType* const image, vtkImageData* const outputImage);

/** Create a color VTK image from the first 3 channels of a vector image. */
void ITKImageToVTKRGBImage(const ImageType* const image, vtkImageData* const outputImage);

/** Create a grayscale VTK image from the magnitude image of a vector image. */
void ITKImageToVTKMagnitudeImage(const ImageType* const image, vtkImageData* const outputImage);

/** Create a grayscale VTK image from the specified 'channel' of the vector image. */
void ITKImageChannelToVTKImage(const ImageType* const image, const unsigned int channel,
                               vtkImageData* const outputImage);

/** Convert a scalar image to a VTK image. */
void ITKScalarImageToVTKImage(const MaskImageType* const image, vtkImageData* const outputImage);

/** Dilate the pixels in 'pixelList' with a kernel of radius 'radius'. */
std::vector<itk::Index<2> > DilatePixelList(const std::vector<itk::Index<2> >& pixelList,
                                            const itk::ImageRegion<2>& region, const unsigned int radius);

/** Convert the points in a polydata to a list of indices. */
std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* const polydata);

void ReplaceChannel(const itk::VectorImage<float, 2>* const image, const unsigned int channel,
                    const itk::Image<float, 2>* const replacement,
                    itk::VectorImage<float, 2>* const output);

/** Compute the median value of a vector. */
template<typename T>
T VectorMedian(const std::vector<T> &v);

/** Compute the average value in a vector. */
template<typename T>
T VectorAverage(const std::vector<T> &v);

/** Write 'image' to 'fileName'. */
template<typename TImage>
void WriteImage(const TImage* const image, const std::string& fileName);

/** Set all pixels in 'pixels' to 'value' in 'image'. */
template <typename TImage>
void SetPixels(TImage* image, const std::vector<itk::Index<2> >& pixels, typename TImage::PixelType& value);

/** Set all pixels in 'pixels' to 'value' in 'image'. */
template <typename TImage>
void SetPixelsInRegionToValue(TImage* const image, const itk::ImageRegion<2>& region,
                              const typename TImage::PixelType& value);

/** Copy an image. */
template<typename TImage>
void DeepCopy(const TImage* const input, TImage* const output);

/** An overload to copy a vector image - the pixel size must be set. */
template<typename TPixel>
void DeepCopy(const itk::VectorImage<TPixel, 2>* const input,
              itk::VectorImage<TPixel, 2>* const output);

/** Extract a channel of an image. */
template<typename TPixel>
void ExtractChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                    typename itk::Image<TPixel, 2>* const output);

template<typename TImage>
void ITKScalarImageToVTKImage(const TImage* const image, vtkImageData* const outputImage);

template<typename TImage>
void NormalizeImage(const TImage* const image, TImage* const outputImage);

template<typename TImage>
float MeanValue(const TImage* const image);

// template<typename TPixel>
// void ReplaceChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
//                     typename itk::Image<TPixel, 2>* const replacement,
//                     typename itk::VectorImage<TPixel, 2>* const output);



template<typename TImage>
float Variance(const TImage* const image);

template<typename TImage>
float StandardDeviation(const TImage* const image);

}

#include "Helpers.hpp"

#endif
