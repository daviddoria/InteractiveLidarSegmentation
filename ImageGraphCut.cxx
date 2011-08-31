/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ImageGraphCut.h"

// ITK
#include "itkBilateralImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

// STL
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <numeric> // accumulate

// VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkLine.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>

// Qt
#include <QMessageBox>

// Custom
#include "Helpers.h"

ImageGraphCut::ImageGraphCut()
{
  this->DifferenceFunction = NULL;
  this->AverageDepthDifference = 0.0;
  this->MedianDepthDifference = 0.0;
  this->MinDepthDifference = 0.0;
  this->MaxDepthDifference = 0.0;
  
  this->AverageColorDifference = 0.0;
  this->MedianColorDifference = 0.0;
  this->MinColorDifference = 0.0;
  this->MaxColorDifference = 0.0;
  
  this->Debug = false;
  
  this->IncludeDepthInHistogram = false;
}

ImageType::Pointer ImageGraphCut::GetImage()
{
  return this->Image;
}

void ImageGraphCut::SetImage(ImageType::Pointer image)
{
  this->Image = ImageType::New();
  this->Image->Graft(image);

  // Setup the output (mask) image
  //this->SegmentMask = GrayscaleImageType::New();
  this->SegmentMask = MaskImageType::New();
  this->SegmentMask->SetRegions(this->Image->GetLargestPossibleRegion());
  this->SegmentMask->Allocate();

  // Setup the image to store the node ids
  this->NodeImage = NodeImageType::New();
  this->NodeImage->SetRegions(this->Image->GetLargestPossibleRegion());
  this->NodeImage->Allocate();

  // Default paramters
  this->Lambda = 0.01;
  this->NumberOfHistogramBins = 10; // This value is never used - it is set from the slider

  // Initializations
  this->ForegroundHistogram = NULL;
  this->BackgroundHistogram = NULL;

  
}

ImageType::Pointer ImageGraphCut::GetMaskedOutput()
{
  // Note: If you get a compiler error on this function complaining about NumericTraits in MaskImageFilter,
  // you will need a newer version of ITK. The ability to mask a VectorImage is new.
  
  // Mask the input image with the mask
  //typedef itk::MaskImageFilter< TImage, GrayscaleImageType > MaskFilterType;
  typedef itk::MaskImageFilter< ImageType, MaskImageType > MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  
  typedef itk::VariableLengthVector<double> VariableVectorType;
  VariableVectorType variableLengthVector;
  variableLengthVector.SetSize(this->Image->GetNumberOfComponentsPerPixel());
  variableLengthVector.Fill(0);
  maskFilter->SetOutsideValue(variableLengthVector);
  
  maskFilter->SetInput1(this->Image);
  maskFilter->SetInput2(this->SegmentMask);
  maskFilter->Update();

  return maskFilter->GetOutput();
}

void ImageGraphCut::CutGraph()
{
  //std::cout << "RGBWeight: " << RGBWeight << std::endl;
  
  // Compute max-flow
  this->Graph->maxflow();

  // Setup the values of the output (mask) image
  //GrayscalePixelType sinkPixel;
  //sinkPixel[0] = 0;
  MaskImageType::PixelType sinkPixel = 0;

  //GrayscalePixelType sourcePixel;
  //sourcePixel[0] = 255;
  MaskImageType::PixelType sourcePixel = 255;

  // Iterate over the node image, querying the Kolmorogov graph object for the association of each pixel and storing them as the output mask
  itk::ImageRegionConstIterator<NodeImageType> nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
    {
    if(this->Graph->what_segment(nodeImageIterator.Get()) == GraphType::SOURCE)
      {
      this->SegmentMask->SetPixel(nodeImageIterator.GetIndex(), sourcePixel);
      }
    else if(this->Graph->what_segment(nodeImageIterator.Get()) == GraphType::SINK)
      {
      this->SegmentMask->SetPixel(nodeImageIterator.GetIndex(), sinkPixel);
      }
    ++nodeImageIterator;
    }

  delete this->Graph;
}

void ImageGraphCut::PerformSegmentation()
{
  // This function performs some initializations and then creates and cuts the graph

  // Ensure at least one pixel has been specified for both the foreground and background
  if((this->Sources.size() <= 0) || (this->Sinks.size() <= 0))
    {
    std::cout << "At least one source (foreground) pixel and one sink (background) pixel must be specified!" << std::endl;
    return;
    }

  // Blank the NodeImage
  itk::ImageRegionIterator<NodeImageType> nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
    {
    nodeImageIterator.Set(NULL);
    ++nodeImageIterator;
    }

  // Blank the output image
  //itk::ImageRegionIterator<GrayscaleImageType> segmentMaskImageIterator(this->SegmentMask, this->SegmentMask->GetLargestPossibleRegion());
  itk::ImageRegionIterator<MaskImageType> segmentMaskImageIterator(this->SegmentMask, this->SegmentMask->GetLargestPossibleRegion());
  segmentMaskImageIterator.GoToBegin();

  MaskImageType::PixelType empty = 0;
  //empty[0] = 0;

  while(!segmentMaskImageIterator.IsAtEnd())
    {
    segmentMaskImageIterator.Set(empty);
    ++segmentMaskImageIterator;
    }

  if(this->Debug)
    {
    this->DifferenceFunction->WriteImages();
    }
  this->CreateGraph();
  this->CutGraph();
}

void ImageGraphCut::CreateHistogram(unsigned int numberOfComponents)
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms
  std::cout << "CreateHistogram()" << std::endl;
  
  // Typedefs
  typedef itk::Statistics::ListSample<PixelType> SampleType;
  typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;

  SampleToHistogramFilterType::Pointer foregroundHistogramFilter = SampleToHistogramFilterType::New();
  SampleToHistogramFilterType::Pointer backgroundHistogramFilter = SampleToHistogramFilterType::New();
  SampleType::Pointer foregroundSample = SampleType::New();
  SampleType::Pointer backgroundSample = SampleType::New();
  
  // We want the histogram bins to take values from 0 to 1 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(numberOfComponents);
  HistogramType::MeasurementVectorType binMaximum(numberOfComponents);
  for(unsigned int i = 0; i < numberOfComponents; i++)
    {
    binMinimum[i] = 0;
    binMaximum[i] = 1;
    }

  // Setup the histogram size
  SampleToHistogramFilterType::HistogramSizeType histogramSize(numberOfComponents);
  histogramSize.Fill(this->NumberOfHistogramBins);

  // Create foreground samples and histogram
  foregroundSample->Clear();
  foregroundSample->SetMeasurementVectorSize(numberOfComponents);
  //std::cout << "Measurement vector size: " << this->ForegroundSample->GetMeasurementVectorSize() << std::endl;
  //std::cout << "Pixel size: " << this->Image->GetPixel(this->Sources[0]).GetNumberOfElements() << std::endl;
  
  for(unsigned int i = 0; i < this->Sources.size(); i++)
    {
    if(!this->Image->GetPixel(this->Sources[i])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
      
    itk::VariableLengthVector<float> normalizedPixel;
    PixelType pixel = this->Image->GetPixel(this->Sources[i]);
    normalizedPixel.SetSize(numberOfComponents);
    for(unsigned int component = 0; component < numberOfComponents; component++)
      {
      normalizedPixel[component] = (pixel[component] - this->DifferenceFunction->MinimumOfChannels[component])/(this->DifferenceFunction->MaximumOfChannels[component] - this->DifferenceFunction->MinimumOfChannels[component]);
      }
    
    foregroundSample->PushBack(normalizedPixel);
    }

  foregroundHistogramFilter->SetHistogramSize(histogramSize);
  foregroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  foregroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  foregroundHistogramFilter->SetAutoMinimumMaximum(false);
  foregroundHistogramFilter->SetInput(foregroundSample);
  foregroundHistogramFilter->Modified();
  foregroundHistogramFilter->Update();
  foregroundHistogramFilter->Register();

  this->ForegroundHistogram = foregroundHistogramFilter->GetOutput();
  //this->ForegroundHistogram->

  // Create background samples and histogram
  backgroundSample->Clear();
  backgroundSample->SetMeasurementVectorSize(numberOfComponents);
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
    {
    if(!this->Image->GetPixel(this->Sinks[i])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
    itk::VariableLengthVector<float> normalizedPixel;
    PixelType pixel = this->Image->GetPixel(this->Sinks[i]);
    normalizedPixel.SetSize(numberOfComponents);
    for(unsigned int component = 0; component < numberOfComponents; component++)
      {
      normalizedPixel[component] = (pixel[component] - this->DifferenceFunction->MinimumOfChannels[component])/(this->DifferenceFunction->MaximumOfChannels[component] - this->DifferenceFunction->MinimumOfChannels[component]);
      }
    backgroundSample->PushBack(normalizedPixel);
    }

  backgroundHistogramFilter->SetHistogramSize(histogramSize);
  backgroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  backgroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  backgroundHistogramFilter->SetAutoMinimumMaximum(false);
  backgroundHistogramFilter->SetInput(backgroundSample);
  backgroundHistogramFilter->Modified();
  backgroundHistogramFilter->Update();
  backgroundHistogramFilter->Register();

  this->BackgroundHistogram = backgroundHistogramFilter->GetOutput();

}


void ImageGraphCut::CreateGraphNodes()
{
  
  // Form the graph
  this->Graph = new GraphType;

  // Add all of the nodes to the graph and store their IDs in a "node image"
  itk::ImageRegionIterator<NodeImageType> nodeImageIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  nodeImageIterator.GoToBegin();

  while(!nodeImageIterator.IsAtEnd())
    {
    nodeImageIterator.Set(this->Graph->add_node());
    ++nodeImageIterator;
    }
}

void ImageGraphCut::CreateNWeights()
{
  
  ////////// Create n-edges and set n-edge weights (links between image nodes) //////////
  // We use a neighborhood iterator here even though we are looking only at a single pixel index in all images on each iteration because we use the neighborhood to determine edge validity.
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  // Traverse the image adding an edge between:
  // - the current pixel and the pixel below it
  // - the current pixel and the pixel to the right of it
  // - the current pixel and the pixel to the bottom-right of it
  // This prevents duplicate edges (i.e. we cannot add an edge to all 8-connected neighbors of every pixel or almost every edge would be duplicated.
  std::cout << "Setting N-Weights..." << std::endl;
  
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();
  
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      //float weight = std::numeric_limits<float>::max(); // This will be the assigned weight if the edge is not computed (if one or both of the pixels is invalid)
      float weight = 0.0;
      bool inbounds = false;
      ImageType::PixelType neighborPixel = iterator.GetPixel(neighbors[i], inbounds);

      // If the current neighbor is outside the image, skip it
      if(!inbounds)
	{
	continue;
	}
	
      // If pixel or its neighbor is not valid, skip this edge.
      if(neighborPixel[4] && centerPixel[4]) // validity channel
	{
	/*
	float depthDifference = depthGradientMagnitudeImage->GetPixel(iterator.GetIndex());
	//float colorDifference = rgbGradientMagnitudeImage->GetPixel(iterator.GetIndex() + neighbors[i]);
	float colorDifference = rgbGradientMagnitudeImage->GetPixel(iterator.GetIndex());
	//float pixelDifference = std::max(depthDifference, colorDifference);
	float pixelDifference = std::max(depthDifference, colorDifference) + (depthDifference + colorDifference)/2.0;
	*/
	
	float pixelDifference = this->DifferenceFunction->GetDifference(iterator.GetIndex());
	  
	// Compute the edge weight
	weight = ComputeEdgeWeight(pixelDifference);

	}// end if current and neighbor are valid
      // Add the edge to the graph
      void* node1 = this->NodeImage->GetPixel(iterator.GetIndex());
      void* node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
      this->Graph->add_edge(node1, node2, weight, weight); // This is an undirected graph so we create a bidirectional edge with both weights set to 'weight'
      //std::cout << "Set n-edge weight to " << weight << std::endl;
      } // end loop over neighbors
    } // end iteration over entire image

}

void ImageGraphCut::CreateTWeights()
{
  
  ////////// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node) //////////

  // Compute the histograms of the selected foreground and background pixels
  //
  //CreateFullHistogramSamples();
  unsigned int numberOfHistogramComponents = 0;
  if(this->IncludeDepthInHistogram)
    {
    numberOfHistogramComponents = 4;
    }
  else
    {
    numberOfHistogramComponents = 3;
    }
    
  CreateHistogram(numberOfHistogramComponents);
  
  itk::ImageRegionIterator<ImageType> imageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NodeImageType> nodeIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  nodeIterator.GoToBegin();

  // Since the t-weight function takes the log of the histogram value,
  // we must handle bins with frequency = 0 specially (because log(0) = -inf)
  // For empty histogram bins we use tinyValue instead of 0.
  float tinyValue = 1e-10;

  std::cout << "Setting T-Weights..." << std::endl;
  
  // These are only for debuging/tracking
  std::vector<float> sinkTWeights;
  std::vector<float> sourceTWeights;
  std::vector<float> sourceHistogramValues;
  std::vector<float> sinkHistogramValues;
  
  // Use the colors only for the t-weights
  while(!imageIterator.IsAtEnd())
    {
    PixelType pixel = imageIterator.Get();
    //float sinkHistogramValue = 0.0;
    //float sourceHistogramValue = 0.0;
    float sinkHistogramValue = tinyValue;
    float sourceHistogramValue = tinyValue;
      
    if(pixel[4])
      {
      //std::cout << "Pixels have size: " << pixel.Size() << std::endl;
      
      HistogramType::MeasurementVectorType measurementVector(numberOfHistogramComponents);
      for(unsigned int i = 0; i < numberOfHistogramComponents; i++)
	{
	//measurementVector[i] = pixel[i];
	measurementVector[i] = (pixel[i] - this->DifferenceFunction->MinimumOfChannels[i])/(this->DifferenceFunction->MaximumOfChannels[i] - this->DifferenceFunction->MinimumOfChannels[i]);
	}

      sinkHistogramValue = this->BackgroundHistogram->GetFrequency(this->BackgroundHistogram->GetIndex(measurementVector));
      sourceHistogramValue = this->ForegroundHistogram->GetFrequency(this->ForegroundHistogram->GetIndex(measurementVector));

      // Convert the histogram value/frequency to make it as if it came from a normalized histogram
      sinkHistogramValue /= this->BackgroundHistogram->GetTotalFrequency();
      sourceHistogramValue /= this->ForegroundHistogram->GetTotalFrequency();

      if(sinkHistogramValue <= 0)
	{
	sinkHistogramValue = tinyValue;
	}
      if(sourceHistogramValue <= 0)
	{
	sourceHistogramValue = tinyValue;
	}
	
      //std::cout << "Setting background weight to: " << -this->Lambda*log(sinkHistogramValue) << std::endl;
      //std::cout << "Setting foreground weight to: " << -this->Lambda*log(sourceHistogramValue) << std::endl;
      
      sinkHistogramValues.push_back(sinkHistogramValue);
      sourceHistogramValues.push_back(sourceHistogramValue);
      sinkTWeights.push_back(-this->Lambda*log(sinkHistogramValue));
      sourceTWeights.push_back(-this->Lambda*log(sourceHistogramValue));
      
      // Add the edge to the graph and set its weight
      this->Graph->add_tweights(nodeIterator.Get(),
                              -this->Lambda*log(sinkHistogramValue),
                              -this->Lambda*log(sourceHistogramValue)); // log() is the natural log
      }
    else
      {
      this->Graph->add_tweights(nodeIterator.Get(), 0, 0);
      }

    ++imageIterator;
    ++nodeIterator;
    }
    
  std::cout << "Average sinkHistogramValue: " << Helpers::VectorAverage<float>(sinkHistogramValues) << std::endl;
  std::cout << "Average sourceHistogramValue: " << Helpers::VectorAverage<float>(sourceHistogramValues) << std::endl;
  
  std::cout << "Max sinkHistogramValue: " << *(std::max_element(sinkHistogramValues.begin(), sinkHistogramValues.end())) << std::endl;
  std::cout << "Max sourceHistogramValue: " << *(std::max_element(sourceHistogramValues.begin(), sourceHistogramValues.end())) << std::endl;
  
  std::cout << "Average sourceTWeights: " << Helpers::VectorAverage<float>(sourceTWeights) << std::endl;
  std::cout << "Average sinkTWeights: " << Helpers::VectorAverage<float>(sinkTWeights) << std::endl;
  
  std::cout << "Max sourceTWeights: " << *(std::max_element(sourceTWeights.begin(), sourceTWeights.end())) << std::endl;
  std::cout << "Max sinkTWeights: " << *(std::max_element(sinkTWeights.begin(), sinkTWeights.end())) << std::endl;
}

void ImageGraphCut::CreateGraph()
{
  if(this->Debug)
    {
    std::cout << "CreateGraph()" << std::endl;
    }
    
  CreateGraphNodes();

  CreateNWeights();
  
  CreateTWeights();
  
  // Set very high source weights for the pixels which were selected as foreground by the user.
  // The syntax is add_tweights(location, SourceLinkWeight, SinkLinkWeight).
  // We want source pixels to be strongly linked to the imaginary source node (so the link won't be cut) and weakly linked
  // to the sink node (weight 0, so it will certainly be cut if necessary. The reverse holds for sink pixels.
  for(unsigned int i = 0; i < this->Sources.size(); i++)
    {
    //this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sources[i]),this->Lambda * std::numeric_limits<float>::max(),0);
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sources[i]), std::numeric_limits<float>::max(),0);
    }

  // Set very high sink weights for the pixels which were selected as background by the user
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
    {
    //this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sinks[i]),0,this->Lambda * std::numeric_limits<float>::max());
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sinks[i]),0, std::numeric_limits<float>::max());
    }
}


std::vector<itk::Index<2> > ImageGraphCut::GetSources()
{
  return this->Sources;
}

void ImageGraphCut::SetLambda(float lambda)
{
  this->Lambda = lambda;
}

void ImageGraphCut::SetNumberOfHistogramBins(int bins)
{
  this->NumberOfHistogramBins = bins;
}

MaskImageType::Pointer ImageGraphCut::GetSegmentMask()
{
  return this->SegmentMask;
}

std::vector<itk::Index<2> > ImageGraphCut::GetSinks()
{
  return this->Sinks;
}

void ImageGraphCut::SetSources(vtkPolyData* sources)
{
  // Convert the vtkPolyData produced by the vtkImageTracerWidget to a list of pixel indices

  this->Sources.clear();

  for(vtkIdType i = 0; i < sources->GetNumberOfPoints(); i++)
    {
    itk::Index<2> index;
    double p[3];
    sources->GetPoint(i,p);
    /*
    index[0] = round(p[0]);
    index[1] = round(p[1]);
    */
    index[0] = vtkMath::Round(p[0]);
    index[1] = vtkMath::Round(p[1]);

    this->Sources.push_back(index);
    }

}

void ImageGraphCut::SetSinks(vtkPolyData* sinks)
{
  // Convert the vtkPolyData produced by the vtkImageTracerWidget to a list of pixel indices

  this->Sinks.clear();

  for(vtkIdType i = 0; i < sinks->GetNumberOfPoints(); i++)
    {
    itk::Index<2> index;
    double p[3];
    sinks->GetPoint(i,p);
    /*
    index[0] = round(p[0]);
    index[1] = round(p[1]);
    */
    index[0] = vtkMath::Round(p[0]);
    index[1] = vtkMath::Round(p[1]);

    this->Sinks.push_back(index);
    }

}

void ImageGraphCut::SetSources(std::vector<itk::Index<2> > sources)
{
  this->Sources = sources;
}

void ImageGraphCut::SetSinks(std::vector<itk::Index<2> > sinks)
{
  this->Sinks = sinks;
}


itk::Size<2> ImageGraphCut::Get1x1Radius()
{
  itk::Size<2> radius;
  radius.Fill(1);
  return radius;
}

void ImageGraphCut::ConstructNeighborhoodIterator(NeighborhoodIteratorType* iterator, std::vector<NeighborhoodIteratorType::OffsetType>& neighbors)
{
  // We are using an 8-connected structure, so the kernel (iteration neighborhood) must only be 3x3 (specified by a radius of 1)
  

  // Traverse the image comparing:
  // - the current pixel and the pixel below it
  // - the current pixel and the pixel to the right of it
  // - the current pixel and the pixel to the bottom-right of it
  // - the current pixel and the pixel to the top-right of it
  
  NeighborhoodIteratorType::OffsetType bottom = {{0,1}};
  neighbors.push_back(bottom);
  NeighborhoodIteratorType::OffsetType right = {{1,0}};
  neighbors.push_back(right);
  NeighborhoodIteratorType::OffsetType bottomRight = {{1,1}};
  neighbors.push_back(bottomRight);
  NeighborhoodIteratorType::OffsetType topRight = {{1,-1}};
  neighbors.push_back(topRight);

  //iterator.Initialize(radius, this->Image, this->Image->GetLargestPossibleRegion());
  
  iterator->ClearActiveList();
  iterator->ActivateOffset(bottom);
  iterator->ActivateOffset(right);
  iterator->ActivateOffset(bottomRight);
  iterator->ActivateOffset(topRight);
}


float ImageGraphCut::ComputeEdgeWeight(float difference)
{
  // Note this value (this->AverageDepthDifference, this->CameraNoise, etc) must correspond to the variance (aka average) of the difference function you are using over the whole image.
  
  
  
  return exp(-pow(difference,2)/(2.0*this->DifferenceFunction->AverageDifference*this->DifferenceFunction->AverageDifference));
}
