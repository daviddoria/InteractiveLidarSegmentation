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
#include "itkImageRegionIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkMaskImageFilter.h"

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

  this->ForegroundSample = SampleType::New();
  this->BackgroundSample = SampleType::New();

  this->ForegroundHistogramFilter = SampleToHistogramFilterType::New();
  this->BackgroundHistogramFilter = SampleToHistogramFilterType::New();
  
  this->CameraNoise = 0.0;
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
  
  // Estimate the "camera noise"
  this->CameraNoise = this->ComputeNoise();

  this->ComputeGlobalStatistics();
  //this->WriteEdges("Edges.vtp");
  this->CreateGraph();
  this->CutGraph();
}

void ImageGraphCut::CreateSamples()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms

  // We want the histogram bins to take values from 0 to 255 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(this->Image->GetNumberOfComponentsPerPixel());
  HistogramType::MeasurementVectorType binMaximum(this->Image->GetNumberOfComponentsPerPixel());
  for(unsigned int i = 0; i < this->Image->GetNumberOfComponentsPerPixel(); i++)
    {
    binMinimum[i] = 0;
    binMaximum[i] = 255;
    }

  // Setup the histogram size
  //std::cout << "Image components per pixel: " << this->Image->GetNumberOfComponentsPerPixel() << std::endl;
  SampleToHistogramFilterType::HistogramSizeType histogramSize(this->Image->GetNumberOfComponentsPerPixel());
  histogramSize.Fill(this->NumberOfHistogramBins);

  // Create foreground samples and histogram
  this->ForegroundSample->Clear();
  this->ForegroundSample->SetMeasurementVectorSize(this->Image->GetNumberOfComponentsPerPixel());
  //std::cout << "Measurement vector size: " << this->ForegroundSample->GetMeasurementVectorSize() << std::endl;
  //std::cout << "Pixel size: " << this->Image->GetPixel(this->Sources[0]).GetNumberOfElements() << std::endl;
  
  for(unsigned int i = 0; i < this->Sources.size(); i++)
    {
    this->ForegroundSample->PushBack(this->Image->GetPixel(this->Sources[i]));
    }

  this->ForegroundHistogramFilter->SetHistogramSize(histogramSize);
  this->ForegroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  this->ForegroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  this->ForegroundHistogramFilter->SetAutoMinimumMaximum(false);
  this->ForegroundHistogramFilter->SetInput(this->ForegroundSample);
  this->ForegroundHistogramFilter->Modified();
  this->ForegroundHistogramFilter->Update();

  this->ForegroundHistogram = this->ForegroundHistogramFilter->GetOutput();

  // Create background samples and histogram
  this->BackgroundSample->Clear();
  this->BackgroundSample->SetMeasurementVectorSize(this->Image->GetNumberOfComponentsPerPixel());
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
    {
    this->BackgroundSample->PushBack(this->Image->GetPixel(this->Sinks[i]));
    }

  this->BackgroundHistogramFilter->SetHistogramSize(histogramSize);
  this->BackgroundHistogramFilter->SetHistogramBinMinimum(binMinimum);
  this->BackgroundHistogramFilter->SetHistogramBinMaximum(binMaximum);
  this->BackgroundHistogramFilter->SetAutoMinimumMaximum(false);
  this->BackgroundHistogramFilter->SetInput(this->BackgroundSample);
  this->BackgroundHistogramFilter->Modified();
  this->BackgroundHistogramFilter->Update();

  this->BackgroundHistogram = BackgroundHistogramFilter->GetOutput();

}

void ImageGraphCut::ComputeGlobalStatistics()
{
  ComputeAllColorDifferences();
  
  // Find min and max
  std::sort(this->AllColorDifferences.begin(), this->AllColorDifferences.end());
  this->MinColorDifference = this->AllColorDifferences[0];
  std::cout << "MinColorDifference: " << this->MinColorDifference << std::endl;
  this->MaxColorDifference = this->AllColorDifferences[this->AllColorDifferences.size()-1];
  std::cout << "MaxColorDifference: " << this->MaxColorDifference << std::endl;
  
  ComputeAllDepthDifferences();
  
  // Find min and max
  std::sort(this->AllDepthDifferences.begin(), this->AllDepthDifferences.end());
  this->MinDepthDifference = this->AllDepthDifferences[0];
  std::cout << "MinDepthDifference: " << this->MinDepthDifference << std::endl;
  this->MaxDepthDifference = this->AllDepthDifferences[this->AllDepthDifferences.size()-1];
  std::cout << "MaxDepthDifference: " << this->MaxDepthDifference << std::endl;
  
  // Compute averages
  this->AverageColorDifference = ComputeAverageColorDifference();
  std::cout << "averageColorDifference: " << this->AverageColorDifference << std::endl;

  this->AverageDepthDifference = ComputeAverageDepthDifference();
  std::cout << "averageDepthDifference: " << this->AverageDepthDifference << std::endl;

  // Compute medians
  this->MedianColorDifference = ComputeMedianColorDifference();
  std::cout << "medianColorDifference: " << this->MedianColorDifference << std::endl;

  this->MedianDepthDifference = ComputeMedianDepthDifference();
  std::cout << "medianDepthDifference: " << this->MedianDepthDifference << std::endl;
    
  this->DifferenceFunction->AverageColorDifference = AverageColorDifference;
  this->DifferenceFunction->MinColorDifference = MinColorDifference;
  this->DifferenceFunction->MaxColorDifference = MaxColorDifference;
  this->DifferenceFunction->MedianColorDifference = MedianColorDifference;
  
  this->DifferenceFunction->AverageDepthDifference = AverageDepthDifference;
  this->DifferenceFunction->MinDepthDifference = MinDepthDifference;
  this->DifferenceFunction->MaxDepthDifference = MaxDepthDifference;
  this->DifferenceFunction->MedianDepthDifference = MedianDepthDifference;
  
}

void ImageGraphCut::CreateGraph()
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
  
  ////////// Create n-edges and set n-edge weights (links between image nodes) //////////
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  // Traverse the image adding an edge between:
  // - the current pixel and the pixel below it
  // - the current pixel and the pixel to the right of it
  // - the current pixel and the pixel to the bottom-right of it
  // This prevents duplicate edges (i.e. we cannot add an edge to all 8-connected neighbors of every pixel or almost every edge would be duplicated.

  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();
  
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      //float weight = std::numeric_limits<float>::max(); // This will be the assigned weight if the edge is not computed (if one or both of the pixels is invalid)
      float weight = 0.0;
      bool inbounds;
      ImageType::PixelType neighborPixel = iterator.GetPixel(neighbors[i], inbounds);

      // If the current neighbor is outside the image, skip it
      if(!inbounds)
	{
	continue;
	}
	
      // If pixel or its neighbor is not valid, skip this edge.
      if(neighborPixel[4] && centerPixel[4]) // validity channel
	{
	//PixelType neighborPixel = iterator.GetPixel(neighbors[i]);

	float pixelDifference = this->DifferenceFunction->Compute(centerPixel, neighborPixel);
	
	if(pixelDifference < 0)
	  {
	  std::cerr << "pixelDifference = " << pixelDifference << " but cannot be negative!" << std::endl;
	  exit(-1);
	  }
	  
	// Compute the edge weight
	weight = ComputeEdgeWeight(pixelDifference);
	//assert(weight >= 0);
	/*
	if(weight < 0)
	  {
	  std::cerr << "pixelDifference = " << pixelDifference << " so then weight = " << weight << std::endl;
	  exit(-1);
	  }
	*/
	}// end if current and neighbor are valid
      // Add the edge to the graph
      void* node1 = this->NodeImage->GetPixel(iterator.GetIndex());
      void* node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
      this->Graph->add_edge(node1, node2, weight, weight); // This is an undirected graph so we create a bidirectional edge with both weights set to 'weight'
      } // end loop over neighbors
    } // end iteration over entire image

  ////////// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node) //////////

  // Compute the histograms of the selected foreground and background pixels
  CreateSamples();

  itk::ImageRegionIterator<ImageType> imageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NodeImageType> nodeIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  nodeIterator.GoToBegin();

  // Since the t-weight function takes the log of the histogram value,
  // we must handle bins with frequency = 0 specially (because log(0) = -inf)
  // For empty histogram bins we use tinyValue instead of 0.
  float tinyValue = 1e-10;

  // Use the colors only for the t-weights
  while(!imageIterator.IsAtEnd())
    {
    PixelType pixel = imageIterator.Get();
    //float sinkHistogramValue = 0.0;
    //float sourceHistogramValue = 0.0;
    float sinkHistogramValue = tinyValue;
    float sourceHistogramValue = tinyValue;
      
    //if(pixel[4])
    if(1)
      {
      //std::cout << "Pixels have size: " << pixel.Size() << std::endl;
      unsigned int measurementVectorSize = std::min(pixel.Size(), 3u);
      
      HistogramType::MeasurementVectorType measurementVector(measurementVectorSize);
      for(unsigned int i = 0; i < measurementVectorSize; i++)
	{
	measurementVector[i] = pixel[i];
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
      }
    else
      {
      }

    // Add the edge to the graph and set its weight
    this->Graph->add_tweights(nodeIterator.Get(),
                              -this->Lambda*log(sinkHistogramValue),
                              -this->Lambda*log(sourceHistogramValue)); // log() is the natural log
    ++imageIterator;
    ++nodeIterator;
    }

  // Set very high source weights for the pixels which were selected as foreground by the user.
  // The syntax is add_tweights(location, SourceLinkWeight, SinkLinkWeight).
  // We want source pixels to be strongly linked to the imaginary source node (so the link won't be cut) and weakly linked to the sink node (weight 0, so it will certainly be cut if necessary. The reverse holds for sink pixels.
  for(unsigned int i = 0; i < this->Sources.size(); i++)
    {
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sources[i]),this->Lambda * std::numeric_limits<float>::max(),0);
    }

  // Set very high sink weights for the pixels which were selected as background by the user
  for(unsigned int i = 0; i < this->Sinks.size(); i++)
    {
    this->Graph->add_tweights(this->NodeImage->GetPixel(this->Sinks[i]),0,this->Lambda * std::numeric_limits<float>::max());
    }
}

/*
float ImageGraphCut::PixelDifference(PixelType a, PixelType b)
{
  // Compute the Euclidean distance between N dimensional pixels
  float difference = 0;

  if(this->Image->GetNumberOfComponentsPerPixel() > 3)
    {
    for(unsigned int i = 0; i < 3; i++)
      {
      difference += (this->RGBWeight / 3.) * pow(a[i] - b[i],2);
      }
    for(unsigned int i = 3; i < this->Image->GetNumberOfComponentsPerPixel(); i++)
      {
      difference += (1 - this->RGBWeight) / (this->Image->GetNumberOfComponentsPerPixel() - 3.) * pow(a[i] - b[i],2);
      }
    }
  else // image is RGB or less (grayscale)
    {
    for(unsigned int i = 0; i < this->Image->GetNumberOfComponentsPerPixel(); i++)
      {
      difference += pow(a[i] - b[i],2);
      }
    }
  return sqrt(difference);
}
*/


/*
float ImageGraphCut::PixelDifferenceUseColorAtSmallDepths(PixelType a, PixelType b)
{
  //float colorDifference = AbsoluteNormalizedPixelColorDifference(a,b);
  float colorDifference = DataNormalizedPixelColorDifference(a,b);
  float depthDifference = DataNormalizedPixelDepthDifference(a,b); // using normalized values is much better because then the main graph cut lambda does not change wildly from data set to data set
  //float depthDifference = PixelDepthDifference(a,b);

  if(depthDifference > 0.4)
    {
    return depthDifference;
    }
  else
    {
    //return depthDifference * (1.0 / colorDifference);
    return depthDifference *  colorDifference;
    }
    return depthDifference;
    
  return depthDifference *  colorDifference;
}
*/


double ImageGraphCut::ComputeNoise()
{
  // Compute an estimate of the "camera noise". This is used in the N-weight function.

  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  double sigma = 0.0;
  int numberOfEdges = 0;

  // Traverse the image collecting the differences between neighboring pixel intensities
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();
    if(!centerPixel[4])
      {
      continue;
      }
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inbounds;
      iterator.GetPixel(neighbors[i], inbounds);
      if(!inbounds)
        {
        continue;
        }
      PixelType neighborPixel = iterator.GetPixel(neighbors[i]);
      if(!neighborPixel[4])
	{
	continue;
	}
      //DifferenceColor differenceColor(*(this->DifferenceFunction));
      DifferenceColorDataNormalized differenceColor(*(this->DifferenceFunction));
      float colorDifference = differenceColor.Compute(centerPixel, neighborPixel);
      sigma += colorDifference;
      numberOfEdges++;
      }

    }

  // Normalize
  sigma /= static_cast<double>(numberOfEdges);

  return sigma;
}


void ImageGraphCut::SetRGBWeight(float weight)
{
  this->RGBWeight = weight;
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

bool ImageGraphCut::IsNaN(const double a)
{
  return a != a;
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

void ImageGraphCut::ComputeAllColorDifferences()
{
  this->AllColorDifferences.clear();
  
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();
    // If pixel is not valid, skip this edge.
    if(!centerPixel[4]) // validity channel
      {
      continue;
      }
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inbounds;
      ImageType::PixelType neighborPixel = iterator.GetPixel(neighbors[i], inbounds);
      
      // If the current neighbor is outside the image, skip it
      if(!inbounds)
        {
        continue;
        }
        
      // If neighbor is not valid, skip this edge.
      if(!neighborPixel[4]) // validity channel
	{
	continue;
	}

      DifferenceColor differenceColor;
      float pixelColorDifference = differenceColor.Compute(centerPixel, neighborPixel);
      this->AllColorDifferences.push_back(pixelColorDifference);
      }
    }

}

float ImageGraphCut::ComputeMedianColorDifference()
{
  float medianColorDifference = Helpers::VectorMedian<float>(this->AllColorDifferences);

  return medianColorDifference;
}

void ImageGraphCut::ComputeAllDepthDifferences()
{
  this->AllDepthDifferences.clear();
  
    std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  std::vector<float> differences;

  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();
    // If pixel is not valid, skip this edge.
    if(!centerPixel[4]) // validity channel
      {
      continue;
      }
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inbounds;
      ImageType::PixelType neighborPixel = iterator.GetPixel(neighbors[i], inbounds);

      // If the current neighbor is outside the image, skip it
      if(!inbounds)
        {
        continue;
        }
      // If pixel is not valid, skip this edge.
      if(!neighborPixel[4]) // validity channel
	{
	continue;
	}

      DifferenceDepth differenceDepth;
      float pixelDifference = differenceDepth.Compute(centerPixel, neighborPixel);

      //std::cout << "pixelDepthDifference: " << pixelDifference << std::endl;
      this->AllDepthDifferences.push_back(pixelDifference);
      }
    }

}

float ImageGraphCut::ComputeMedianDepthDifference()
{
  float medianDifference = Helpers::VectorMedian<float>(this->AllDepthDifferences);

  return medianDifference;
}

float ImageGraphCut::ComputeAverageColorDifference()
{
  float sumColorDifferences = std::accumulate(this->AllColorDifferences.begin(), this->AllColorDifferences.end(), 0);
  float averageColorDifference = sumColorDifferences / static_cast<float>(this->AllColorDifferences.size());

  return averageColorDifference;
}

float ImageGraphCut::ComputeAverageDepthDifference()
{
  float sumDepthDifferences = std::accumulate(this->AllDepthDifferences.begin(), this->AllDepthDifferences.end(), 0);
  float averageDepthDifference = sumDepthDifferences / static_cast<float>(this->AllDepthDifferences.size());

  return averageDepthDifference;
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


void ImageGraphCut::WriteEdges(const std::string& fileName)
{
  UsedDepth = 0;
  UsedColor = 0;
  
  // Create a vtkPoints object that we will later create line segments on. Store the PointIDs in an image for easy lookup (corresponding to pixels that will be traversed in other images)
  typedef itk::Image<unsigned int, 2> UnsignedIntImageType;
  UnsignedIntImageType::Pointer pointIDImage = UnsignedIntImageType::New();
  pointIDImage->SetRegions(this->Image->GetLargestPossibleRegion());
  pointIDImage->Allocate();
  pointIDImage->FillBuffer(0);
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  
  itk::ImageRegionIterator<UnsignedIntImageType> pointIDImageIterator(pointIDImage, pointIDImage->GetLargestPossibleRegion());
 
  unsigned int counter = 0;
  while(!pointIDImageIterator.IsAtEnd())
    {
    double p[3];
    p[0] = pointIDImageIterator.GetIndex()[0];
    p[1] = pointIDImageIterator.GetIndex()[1];
    p[2] = 0;
    points->InsertNextPoint(p);
    
    pointIDImageIterator.Set(counter);
  
    counter++;
    ++pointIDImageIterator;
    }
  
  
  ////////// Create n-edges and set n-edge weights (links between image nodes) //////////
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
  ConstructNeighborhoodIterator(&iterator, neighbors);

  // Traverse the image adding an edge between:
  // - the current pixel and the pixel below it
  // - the current pixel and the pixel to the right of it
  // - the current pixel and the pixel to the bottom-right of it
  // This prevents duplicate edges (i.e. we cannot add an edge to all 8-connected neighbors of every pixel or almost every edge would be duplicated.
  
  // Create a cell array to store the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  
  vtkSmartPointer<vtkFloatArray> finalWeights = vtkSmartPointer<vtkFloatArray>::New();
  finalWeights->SetNumberOfComponents(1);
  finalWeights->SetName("FinalWeights");
  
  // Color
  vtkSmartPointer<vtkFloatArray> colorWeights = vtkSmartPointer<vtkFloatArray>::New();
  colorWeights->SetNumberOfComponents(1);
  colorWeights->SetName("ColorWeights");
  
  vtkSmartPointer<vtkFloatArray> absoluteNormalizedColorWeights = vtkSmartPointer<vtkFloatArray>::New();
  absoluteNormalizedColorWeights->SetNumberOfComponents(1);
  absoluteNormalizedColorWeights->SetName("AbsoluteNormalizedColorWeights");
  
  vtkSmartPointer<vtkFloatArray> dataNormalizedColorWeights = vtkSmartPointer<vtkFloatArray>::New();
  dataNormalizedColorWeights->SetNumberOfComponents(1);
  dataNormalizedColorWeights->SetName("DataNormalizedColorWeights");
  
  // Depth
  vtkSmartPointer<vtkFloatArray> depthWeights = vtkSmartPointer<vtkFloatArray>::New();
  depthWeights->SetNumberOfComponents(1);
  depthWeights->SetName("DepthWeights");
  
  vtkSmartPointer<vtkFloatArray> absoluteNormalizedDepthWeights = vtkSmartPointer<vtkFloatArray>::New();
  absoluteNormalizedDepthWeights->SetNumberOfComponents(1);
  absoluteNormalizedDepthWeights->SetName("AbsoluteNormalizedDepthWeights");

  vtkSmartPointer<vtkFloatArray> dataNormalizedDepthWeights = vtkSmartPointer<vtkFloatArray>::New();
  dataNormalizedDepthWeights->SetNumberOfComponents(1);
  dataNormalizedDepthWeights->SetName("DataNormalizedDepthWeights");
  
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
    PixelType centerPixel = iterator.GetCenterPixel();

    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inbounds;
      iterator.GetPixel(neighbors[i], inbounds);

      // If the current neighbor is outside the image, skip it
      if(!inbounds)
        {
        continue;
        }
      
      unsigned int currentPixelId = pointIDImage->GetPixel(iterator.GetIndex());
      
      unsigned int neighborPixelId = pointIDImage->GetPixel(iterator.GetIndex(neighbors[i]));
      
      PixelType neighborPixel = iterator.GetPixel(neighbors[i]);
      
      if(centerPixel[4] && neighborPixel[4])
	{

	// Compute the Euclidean distance between the pixel intensities
	DifferenceDepthWeightedByColor differenceDepthWeightedByColor(*(this->DifferenceFunction));
	float finalDifference = differenceDepthWeightedByColor.Compute(centerPixel, neighborPixel);
	
	finalWeights->InsertNextValue(finalDifference);
	
	// Color
	DifferenceColor colorDifferenceFunction(*(this->DifferenceFunction));
	
	float colorDifference = colorDifferenceFunction.Compute(centerPixel, neighborPixel);
	colorWeights->InsertNextValue(colorDifference);
	
	DifferenceColorAbsoluteNormalized absoluteNormalizedPixelColorDifferenceFunction(*(this->DifferenceFunction));
	float absoluteNormalizedColorDifference = absoluteNormalizedPixelColorDifferenceFunction.Compute(centerPixel, neighborPixel);
	absoluteNormalizedColorWeights->InsertNextValue(absoluteNormalizedColorDifference);
	
	DifferenceColorDataNormalized dataNormalizedPixelColorDifferenceFunction(*(this->DifferenceFunction));
	float dataNormalizedColorDifference = dataNormalizedPixelColorDifferenceFunction.Compute(centerPixel, neighborPixel);
	dataNormalizedColorWeights->InsertNextValue(dataNormalizedColorDifference);
	
	// Depth
	DifferenceDepth depthDifferenceFunction(*(this->DifferenceFunction));
	float depthDifference = depthDifferenceFunction.Compute(centerPixel, neighborPixel);
	depthWeights->InsertNextValue(depthDifference);
	
	DifferenceDepthAbsoluteNormalized absoluteNormalizedPixelDepthDifferenceFunction(*(this->DifferenceFunction));
	float absoluteNormalizedDepthDifference = absoluteNormalizedPixelDepthDifferenceFunction.Compute(centerPixel, neighborPixel);
	absoluteNormalizedDepthWeights->InsertNextValue(absoluteNormalizedDepthDifference);
	
	DifferenceDepthDataNormalized dataNormalizedPixelDepthDifferenceFunction(*(this->DifferenceFunction));
	float dataNormalizedDepthDifference = dataNormalizedPixelDepthDifferenceFunction.Compute(centerPixel, neighborPixel);
	dataNormalizedDepthWeights->InsertNextValue(dataNormalizedDepthDifference);
	
	//std::cout << "weight between " << currentPixelId << " and " << neighborPixelId << " set to " << pixelDifference << std::endl;
	}
      else
	{
	finalWeights->InsertNextValue(0);
	colorWeights->InsertNextValue(0);
	absoluteNormalizedColorWeights->InsertNextValue(0);
	dataNormalizedColorWeights->InsertNextValue(0);
	depthWeights->InsertNextValue(0);
	absoluteNormalizedDepthWeights->InsertNextValue(0);
	dataNormalizedDepthWeights->InsertNextValue(0);
	}
      // Create an edge
      vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
      line->GetPointIds()->SetId(0,currentPixelId);
      line->GetPointIds()->SetId(1,neighborPixelId);
      lines->InsertNextCell(line);
 
      } // end neighbors loop
    
    } // end iterator loop

  // Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->GetCellData()->SetScalars(finalWeights);
  polyData->GetCellData()->AddArray(colorWeights);
  polyData->GetCellData()->AddArray(absoluteNormalizedColorWeights);
  polyData->GetCellData()->AddArray(dataNormalizedColorWeights);
  polyData->GetCellData()->AddArray(depthWeights);
  polyData->GetCellData()->AddArray(absoluteNormalizedDepthWeights);
  polyData->GetCellData()->AddArray(dataNormalizedDepthWeights);
  
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(polyData);
  writer->Write();
}

float ImageGraphCut::ComputeEdgeWeight(float difference)
{
  return exp(-pow(difference,2)/(2.0*this->CameraNoise*this->CameraNoise));  
}
