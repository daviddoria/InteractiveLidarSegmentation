#include "Helpers.h"

// ITK
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

// VTK
#include <vtkPolyData.h>

// STL
#include <iostream>

namespace Helpers
{

float NegativeLog(const float value)
{
  return -1. * log(value);
}

std::vector<itk::Index<2> > GetNonZeroPixels(const MaskImageType* const image)
{
  std::vector<itk::Index<2> > pixels;
  itk::ImageRegionConstIterator<MaskImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // If the current pixel is in question
      {
      pixels.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }
  return pixels;
}

unsigned int CountNonZeroPixels(const MaskImageType* const image)
{
  std::vector<itk::Index<2> > pixels = GetNonZeroPixels(image);
  return pixels.size();
}
      
bool FindClosestNonZeroPixel(const MaskImageType* const image, const itk::Index<2>& queryPixel, const unsigned int radiusValue, itk::Index<2>& returnPixel)
{
  itk::Index<2> zeroIndex;
  zeroIndex.Fill(0);
  
  ImageType::SizeType radius;
  radius.Fill(radiusValue);
  
  itk::Size<2> size;
  size.Fill(1);
  itk::ImageRegion<2> region(queryPixel,size);
  
  itk::ConstNeighborhoodIterator<MaskImageType> iterator(radius, image, region);
  
  unsigned int length = (radiusValue*2)+1;
  while(!iterator.IsAtEnd())
    {
    for(unsigned int i = 0; i < length*length; i++)
      {
      if(iterator.GetIndex(i) == queryPixel) // Skip the pixel we are currently at
	{
	continue;
	}
      bool inBounds;
      ImageType::PixelType pixel = iterator.GetPixel(i, inBounds);
      if(pixel != 0) // We found a non-zero pixel
	{
	returnPixel = iterator.GetIndex(i);
	return true;
	}
 
      }
    ++iterator;
    }
    
  return false;
}

itk::Index<2> FindClosestNonZeroPixel(const MaskImageType* const image, const itk::Index<2>& queryPixel)
{
  // Look in successively bigger neighborhoods
  for(unsigned int radiusValue = 1; radiusValue < std::max(image->GetLargestPossibleRegion().GetSize()[0], image->GetLargestPossibleRegion().GetSize()[1]); ++radiusValue)
    {
    //std::cout << "Radius: " << radiusValue << std::endl;
    itk::Index<2> closestPixel;
    bool success = FindClosestNonZeroPixel(image, queryPixel, radiusValue, closestPixel);
    if(success)
      {
      return closestPixel;
      }
    }
  std::cerr << "No non-zero pixel was found!" << std::endl;
  
  itk::Index<2> zeroIndex;
  zeroIndex.Fill(0);
  return zeroIndex;
}
  
bool IsNaN(const double a)
{
  return a != a;
}

void IndicesToBinaryImage(const std::vector<itk::Index<2> >& indices, UnsignedCharScalarImageType* const image)
{
  // The Regions of the 'image' must be set before calling this function
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  image->Allocate();
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    image->SetPixel(indices[i], 255);
    }
}

std::vector<itk::Index<2> > BinaryImageToIndices(const UnsignedCharScalarImageType* const image)
{
  std::vector<itk::Index<2> > indices;
  
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(image,image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      indices.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }

  return indices;
}

void MaskImage(vtkImageData* const VTKImage, vtkImageData* const VTKSegmentMask, vtkImageData* const VTKMaskedImage)
{
  int* dims = VTKImage->GetDimensions();

  VTKMaskedImage->SetDimensions(dims);
  VTKMaskedImage->SetNumberOfScalarComponents(4);
  VTKMaskedImage->SetScalarTypeToUnsignedChar();

  // int dims[3]; // can't do this
  for (int y = 0; y < dims[1]; y++)
    {
    for (int x = 0; x < dims[0]; x++)
      {

      unsigned char* imagePixel = static_cast<unsigned char*>(VTKImage->GetScalarPointer(x,y,0));
      unsigned char* maskPixel = static_cast<unsigned char*>(VTKSegmentMask->GetScalarPointer(x,y,0));
      unsigned char* outputPixel = static_cast<unsigned char*>(VTKMaskedImage->GetScalarPointer(x,y,0));

      outputPixel[0] = imagePixel[0];

      if(VTKImage->GetNumberOfScalarComponents() == 3)
        {
        outputPixel[1] = imagePixel[1];
        outputPixel[2] = imagePixel[2];
        }
      else // Grayscale should have all components equal to the first component
        {
        outputPixel[1] = imagePixel[0];
        outputPixel[2] = imagePixel[0];
        }

      if(maskPixel[0] == 0)
        {
        outputPixel[3] = 0;
        }
      else
        {
        outputPixel[3] = 255;
        }

      }
    }
}


// Convert single channel ITK image to VTK image
void ITKScalarImageToVTKImage(const MaskImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKScalarImagetoVTKImage()" << std::endl;
  
  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(1);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<MaskImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = static_cast<unsigned char>(imageIterator.Get());
    ++imageIterator;
    }
}


// Convert a vector ITK image to a VTK image for display
void ITKImageToVTKImage(const ImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "Enter ITKImagetoVTKImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() >= 3)
    {
    ITKImageToVTKRGBImage(image, outputImage);
    }
  else
    {
    ITKImageToVTKMagnitudeImage(image, outputImage);
    }
  //std::cout << "Exit ITKImagetoVTKImage()" << std::endl;
}

// Convert a vector ITK image to a VTK image for display
void ITKImageToVTKRGBImage(const ImageType* const image, vtkImageData* const outputImage)
{
  // This function assumes an ND (with N>3) image has the first 3 channels as RGB and
  // extra information in the remaining channels.
  
  //std::cout << "Enter ITKImagetoVTKRGBImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() < 3)
    {
    std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel()
              << " components, but at least 3 are required." << std::endl;
    return;
    }

  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(3);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<ImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < 3; component++)
      {
      pixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
      }

    ++imageIterator;
    }
    
  //std::cout << "Exit ITKImagetoVTKRGBImage()" << std::endl;
}


// Convert a vector ITK image to a VTK image for display
void ITKImageToVTKMagnitudeImage(const ImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKImagetoVTKMagnitudeImage()" << std::endl;
  // Compute the magnitude of the ITK image
  typedef itk::VectorMagnitudeImageFilter<
                  ImageType, FloatScalarImageType >  VectorMagnitudeFilterType;

  // Create and setup a magnitude filter
  VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
  magnitudeFilter->SetInput( image );
  magnitudeFilter->Update();

  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<
                  FloatScalarImageType, UnsignedCharScalarImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput( magnitudeFilter->GetOutput() );
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  outputImage->SetNumberOfScalarComponents(1);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType>
    imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());

  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }
}

void ITKImageChannelToVTKImage(const ImageType* const image, const unsigned int channel,
                               vtkImageData* const outputImage)
{
  if(channel >= image->GetNumberOfComponentsPerPixel())
  {
    std::cerr << "Cannot extract channel " << channel << " of a "
              << image->GetNumberOfComponentsPerPixel() << " channel image." << std::endl;
    return;
  }

  typedef itk::Image<typename ImageType::InternalPixelType, 2> ScalarImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<ImageType, ScalarImageType> IndexSelectionType;
  IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(channel);
  indexSelectionFilter->SetInput(image);
  indexSelectionFilter->Update();

  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<
                  ScalarImageType, UnsignedCharScalarImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(indexSelectionFilter->GetOutput());
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  outputImage->SetNumberOfScalarComponents(1);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType>
    imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());

  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }
}

void InvertBinaryImage(const UnsignedCharScalarImageType* const image, UnsignedCharScalarImageType* const inverted)
{
  // Setup inverted image
  inverted->SetRegions(image->GetLargestPossibleRegion());
  inverted->Allocate();
  
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(image,image->GetLargestPossibleRegion());
 
  while(!imageIterator.IsAtEnd())
    {
    // Get the value of the current pixel
    if(imageIterator.Get())
      {
      inverted->SetPixel(imageIterator.GetIndex(), 0);
      }
    else
      {
      inverted->SetPixel(imageIterator.GetIndex(), 255);
      }
 
    ++imageIterator;
    }
}


std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* const polydata)
{
  // The points of the polydata are floating point values, we must convert them to pixel indices.
  
  //std::cout << "Enter PolyDataToPixelList()" << std::endl;
  std::cout << "There are " << polydata->GetNumberOfPoints() << " points." << std::endl;
  
  // Convert vtkPoints to indices
  //std::cout << "Converting vtkPoints to indices..." << std::endl;
  std::vector<itk::Index<2> > linePoints;
  for(vtkIdType pointId = 0; pointId < polydata->GetNumberOfPoints(); ++pointId)
    {
    itk::Index<2> index;
    double p[3];
    polydata->GetPoint(pointId, p);
    // std::cout << "point " << pointId << " : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    index[0] = round(p[0]);
    index[1] = round(p[1]);
    if(linePoints.size() == 0)
      {
      linePoints.push_back(index);
      continue;
      }

    // Don't duplicate indices of points acquired in a row that round to the same pixel.
    if(index != linePoints[linePoints.size() - 1])
      {
      linePoints.push_back(index);
      }
    }

  if(linePoints.size() < 2)
    {
    std::cerr << "Cannot draw a lines between " << linePoints.size() << " points." << std::endl;
    return linePoints;
    }
    
  // Compute the indices between every pair of points
  //std::cout << "Computing the indices between every pair of points..." << std::endl;
  std::vector<itk::Index<2> > allIndices;
  for(unsigned int linePointId = 1; linePointId < linePoints.size(); linePointId++)
    {
    //std::cout << "Getting the indices..." << std::endl;
    itk::Index<2> index0 = linePoints[linePointId-1];
    itk::Index<2> index1 = linePoints[linePointId];

    if(index0 == index1)
      {
      std::cout << "Can't draw a line between the same pixels (" << index0 << " and " << index1 << "!" << std::endl;
      continue;
      }

    //std::cout << "Constructing the line..." << std::endl;
    itk::BresenhamLine<2> line;
    std::vector<itk::Index<2> > indices = line.BuildLine(index0, index1);
    //std::cout << "Saving indices..." << std::endl;
    for(unsigned int i = 0; i < indices.size(); i++)
      {
      allIndices.push_back(indices[i]);
      }
      
    } // end for loop over line segments

  //std::cout << "Exit PolyDataToPixelList()" << std::endl;
  return allIndices;
}

std::vector<itk::Index<2> > DilatePixelList(const std::vector<itk::Index<2> >& pixelList,
                                            const itk::ImageRegion<2>& region, const unsigned int radius)
{
  //std::cout << "DilatePixelList: input has " << pixelList.size() << " pixels." << std::endl;
  // Construct an image of the pixels in the list
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  typedef std::vector<itk::Index<2> > PixelVectorType;
  
  for(PixelVectorType::const_iterator iter = pixelList.begin(); iter != pixelList.end(); ++iter)
  {
    // Note, this must be 255, not just any non-zero number, for BinaryDilateImageFilter to work properly.
    image->SetPixel(*iter, 255); 
  }

  //WriteImage(image.GetPointer(), "beforeDilation.png");

  // Dilate the image
  typedef itk::BinaryBallStructuringElement<ImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(image);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  //WriteImage(dilateFilter->GetOutput(), "afterDilation.png");
  
  PixelVectorType dilatedPixelList;
  
  itk::ImageRegionConstIteratorWithIndex<ImageType> imageIterator(dilateFilter->GetOutput(),
                                                         dilateFilter->GetOutput()->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      dilatedPixelList.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }

  //std::cout << "DilatePixelList: output has " << dilatedPixelList.size() << " pixels." << std::endl;
  return dilatedPixelList;
}

void CreateTransparentImage(vtkImageData* const VTKImage)
{
  VTKImage->SetNumberOfScalarComponents(4);
  VTKImage->SetScalarTypeToUnsignedChar();
  VTKImage->AllocateScalars();
  
  int* dims = VTKImage->GetDimensions();
  // int dims[3]; // can't do this
  for (int y = 0; y < dims[1]; y++)
    {
    for (int x = 0; x < dims[0]; x++)
      {
      unsigned char* outputPixel = static_cast<unsigned char*>(VTKImage->GetScalarPointer(x,y,0));
      unsigned char color = 255;
      outputPixel[0] = color;
      outputPixel[1] = color;
      outputPixel[2] = color;

      outputPixel[3] = 0;
      //outputPixel[3] = 255;
      
      } // end x loop
    } // end y loop
}

void SetPixels(vtkImageData* const VTKImage, const std::vector<itk::Index<2> >& pixels, const unsigned char color[3])
{
  int* dims = VTKImage->GetDimensions();
  
  for(unsigned int i = 0; i < pixels.size(); ++i)
    {
    if(pixels[i][0] >= dims[0] || pixels[i][1] >= dims[1]) // The pixel is out of bounds
      {
      continue;
      }
    unsigned char* pixel = static_cast<unsigned char*>(VTKImage->GetScalarPointer(pixels[i][0],pixels[i][1],0));
    pixel[0] = color[0];
    pixel[1] = color[1];
    pixel[2] = color[2];
    // Make sure the pixel is not transparent
    if(VTKImage->GetNumberOfScalarComponents() == 4)
      {
      pixel[3] = 255;
      }
    }
  
}

void SetImageSize(vtkImageData* input, vtkImageData* output)
{
  int* dims = input->GetDimensions();
  output->SetDimensions(dims); 
}

} // end namespace
