#include "Helpers.h"

// ITK
#include "itkBresenhamLine.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"

// VTK
#include <vtkPolyData.h>

namespace Helpers
{

unsigned int CountNonZeroPixels(MaskImageType::Pointer image)
{
  unsigned int counter = 0;
  itk::ImageRegionIterator<MaskImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // If the current pixel is in question
      {
      counter++;
      }
    ++imageIterator;
    }
  return counter;
}
      
bool FindClosestNonZeroPixel(MaskImageType::Pointer image, itk::Index<2> queryPixel, unsigned int radiusValue, itk::Index<2>& returnPixel)
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

itk::Index<2> FindClosestNonZeroPixel(MaskImageType::Pointer image, itk::Index<2> queryPixel)
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

void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image)
{
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  // Blank the image
  /*
  itk::ImageRegionIterator<UnsignedCharScalarImageType> imageIterator(image,image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(0);
    ++imageIterator;
    }
  */
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    image->SetPixel(indices[i], 255);
    }
}

std::vector<itk::Index<2> > BinaryImageToIndices(UnsignedCharScalarImageType::Pointer image)
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

void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage)
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
void ITKScalarImagetoVTKImage(MaskImageType::Pointer image, vtkImageData* outputImage)
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
void ITKImagetoVTKImage(ImageType::Pointer image, vtkImageData* outputImage)
{
  //std::cout << "Enter ITKImagetoVTKImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() >= 3)
    {
    ITKImagetoVTKRGBImage(image, outputImage);
    }
  else
    {
    ITKImagetoVTKMagnitudeImage(image, outputImage);
    }
  //std::cout << "Exit ITKImagetoVTKImage()" << std::endl;
}

// Convert a vector ITK image to a VTK image for display
void ITKImagetoVTKRGBImage(ImageType::Pointer image, vtkImageData* outputImage)
{
  // This function assumes an ND (with N>3) image has the first 3 channels as RGB and extra information in the remaining channels.
  
  //std::cout << "Enter ITKImagetoVTKRGBImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() < 3)
    {
    std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel() << " components, but at least 3 are required." << std::endl;
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
void ITKImagetoVTKMagnitudeImage(ImageType::Pointer image, vtkImageData* outputImage)
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
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType> imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }
}

void InvertBinaryImage(UnsignedCharScalarImageType::Pointer image, UnsignedCharScalarImageType::Pointer inverted)
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


std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata)
{
  std::cout << "Enter PolyDataToPixelList()" << std::endl;
  std::cout << "There are " << polydata->GetNumberOfPoints() << " points." << std::endl;
  
  // Convert vtkPoints to indices
  std::cout << "Converting vtkPoints to indices..." << std::endl;
  std::vector<itk::Index<2> > linePoints;
  for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
    itk::Index<2> index;
    double p[3];
    polydata->GetPoint(i,p);
    index[0] = round(p[0]);
    index[1] = round(p[1]);
    linePoints.push_back(index);
    }

  if(linePoints.size() < 2)
    {
    std::cerr << "Cannot draw a lines between " << linePoints.size() << " points." << std::endl;
    return linePoints;
    }
    
  // Compute the indices between every pair of points
  std::cout << "Computing the indices between every pair of points..." << std::endl;
  std::vector<itk::Index<2> > allIndices;
  for(unsigned int linePointId = 1; linePointId < linePoints.size(); linePointId++)
    {
    std::cout << "Getting the indices..." << std::endl;
    itk::Index<2> index0 = linePoints[linePointId-1];
    itk::Index<2> index1 = linePoints[linePointId];
    /*
    // Currently need the distance between the points for Bresenham (pending patch in Gerrit)
    itk::Point<float,2> point0;
    itk::Point<float,2> point1;
    for(unsigned int i = 0; i < 2; i++)
      {
      point0[i] = index0[i];
      point1[i] = index1[i];
      }
      
    float distance = point0.EuclideanDistanceTo(point1);
    itk::BresenhamLine<2> line;
    std::vector<itk::Offset<2> > offsets;

    // Occasionally there is a crash when scribbling. One time while attempting to stepp through another part of the code it seemed to crash in BresenhamLine
    // so I added this try/catch to see if it is ever thrown in the future.
    try
    {
      // This indicates to draw a line of length 'distance' starting at (0,0) in the direction (point1-point0). We will add these offsets to index0 later.
      offsets = line.BuildLine(point1-point0, distance);
    }
    catch (...)
    {
      cerr << "Error building BresenhamLine." << endl;
    }
    
    for(unsigned int i = 0; i < offsets.size(); i++)
      {
      allIndices.push_back(index0 + offsets[i]);
      }
    */
    std::cout << "Constructing the line..." << std::endl;
    itk::BresenhamLine<2> line;
    std::vector<itk::Index<2> > indices = line.BuildLine(index0, index1);
    std::cout << "Saving indices..." << std::endl;
    for(unsigned int i = 0; i < indices.size(); i++)
      {
      allIndices.push_back(indices[i]);
      }
      
    } // end for loop over line segments

  std::cout << "Exit PolyDataToPixelList()" << std::endl;
  return allIndices;
}


} // end namespace
