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
//#include "itkAndImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkXorImageFilter.h"

// STL
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <numeric> // for accumulate()

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

  this->Debug = false;
  
  this->IncludeDepthInHistogram = false;
  this->NumberOfHistogramComponents = 0;
  
  // Debug
  this->DebugGraphPolyData = vtkSmartPointer<vtkPolyData>::New();
  this->DebugGraphLines = vtkSmartPointer<vtkCellArray>::New();
  
  this->DebugGraphEdgeWeights = vtkSmartPointer<vtkFloatArray>::New();
  this->DebugGraphEdgeWeights->SetNumberOfComponents(1);
  this->DebugGraphEdgeWeights->SetName("EdgeWeights");
  
  this->DebugGraphSourceWeights = vtkSmartPointer<vtkFloatArray>::New();
  this->DebugGraphSourceWeights->SetNumberOfComponents(1);
  this->DebugGraphSourceWeights->SetName("SourceWeights");
  
  this->DebugGraphSinkWeights = vtkSmartPointer<vtkFloatArray>::New();
  this->DebugGraphSinkWeights->SetNumberOfComponents(1);
  this->DebugGraphSinkWeights->SetName("SinkWeights");
  
  this->DebugGraphSourceHistogram = vtkSmartPointer<vtkFloatArray>::New();
  this->DebugGraphSourceHistogram->SetNumberOfComponents(1);
  this->DebugGraphSourceHistogram->SetName("SourceHistogram");
  
  this->DebugGraphSinkHistogram = vtkSmartPointer<vtkFloatArray>::New();
  this->DebugGraphSinkHistogram->SetNumberOfComponents(1);
  this->DebugGraphSinkHistogram->SetName("SinkHistogram");
  
  this->DebugGraphPointIds = itk::Image<unsigned int, 2>::New();
  
}

void ImageGraphCut::CreateDebugPolyData()
{
  this->DebugGraphPointIds->SetRegions(this->Image->GetLargestPossibleRegion());
  this->DebugGraphPointIds->Allocate();
  this->DebugGraphPointIds->FillBuffer(0);
  
  itk::ImageRegionIterator<itk::Image<unsigned int, 2> > imageIterator(this->DebugGraphPointIds, this->DebugGraphPointIds->GetLargestPossibleRegion());
 
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  
  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(points->GetNumberOfPoints());
    double p[3];
    p[0] = imageIterator.GetIndex()[0];
    p[1] = imageIterator.GetIndex()[1];
    p[2] = 0;
  
    points->InsertNextPoint(p);
    ++imageIterator;
    }
  this->DebugGraphPolyData->SetPoints(points);
}

ImageType::Pointer ImageGraphCut::GetImage()
{
  return this->Image;
}

void ImageGraphCut::SetImage(const ImageType* const image)
{
  this->Image = ImageType::New();
  Helpers::DeepCopy(image, this->Image.GetPointer());

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
  if(this->Debug)
    {
    std::cout << "CutGraph()" << std::endl;
    }
  
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

  // Only keep the largest segment
  typedef itk::ConnectedComponentImageFilter<MaskImageType, MaskImageType> ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New ();
  connectedComponentFilter->SetInput(SegmentMask);
  connectedComponentFilter->Update();

  //std::cout << "Number of objects: " << connectedComponentFilter->GetObjectCount() << std::endl;

  typedef itk::LabelShapeKeepNObjectsImageFilter<MaskImageType> LabelShapeKeepNObjectsImageFilterType;
  LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
  labelShapeKeepNObjectsImageFilter->SetInput(connectedComponentFilter->GetOutput());
  labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
  labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
  labelShapeKeepNObjectsImageFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
  labelShapeKeepNObjectsImageFilter->Update();

  typedef itk::RescaleIntensityImageFilter<MaskImageType, MaskImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
  rescaleFilter->Update();

  Helpers::DeepCopy(rescaleFilter->GetOutput(), SegmentMask.GetPointer());
}

void ImageGraphCut::PerformSegmentation()
{
  std::cout << "PerformSegmentation() " << std::endl;
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

  if(this->IncludeDepthInHistogram)
    {
    this->NumberOfHistogramComponents = 4;
    }
  else
    {
    this->NumberOfHistogramComponents = 3;
    }
    
  if(this->Debug)
    {
    //this->DifferenceFunction->WriteImages();
    }
  this->CreateGraph();
  
  this->CutGraph();
  
  delete this->Graph;
}

const HistogramType* ImageGraphCut::CreateHistogram(std::vector<itk::Index<2> > pixels, std::vector<unsigned int> channelsToUse)
//void ImageGraphCut::CreateHistogram(std::vector<itk::Index<2> > pixels, std::vector<unsigned int> channelsToUse, const HistogramType* histogramOutput)
{
  std::cout << "CreateHistogram()" << std::endl;
  unsigned int numberOfComponents = channelsToUse.size();

  // Typedefs
  typedef itk::Statistics::ListSample<PixelType> SampleType;
  typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;

  SampleToHistogramFilterType::Pointer histogramFilter = SampleToHistogramFilterType::New();

  SampleType::Pointer sample = SampleType::New();
  
  std::vector<float> debugNormalizedPixelValues;
  
  // We want the histogram bins to take values from 0 to 1 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(numberOfComponents);
  HistogramType::MeasurementVectorType binMaximum(numberOfComponents);
  for(unsigned int component = 0; component < numberOfComponents; component++)
    {
    binMinimum[component] = 0;
    binMaximum[component] = 1;
    }

  // Setup the histogram size
  SampleToHistogramFilterType::HistogramSizeType histogramSize(numberOfComponents);
  histogramSize.Fill(this->NumberOfHistogramBins);

  // Create samples and histogram
  sample->Clear();
  sample->SetMeasurementVectorSize(numberOfComponents);
  //std::cout << "Measurement vector size: " << this->ForegroundSample->GetMeasurementVectorSize() << std::endl;
  //std::cout << "Pixel size: " << this->Image->GetPixel(this->Sources[0]).GetNumberOfElements() << std::endl;

  std::vector<float> minimumOfChannels = Helpers::ComputeMinOfAllChannels(this->Image);
  std::vector<float> maximumOfChannels = Helpers::ComputeMaxOfAllChannels(this->Image);
  
  for(unsigned int pixelId = 0; pixelId < pixels.size(); pixelId++) // Add all of the indicated foreground pixels to the histogram
    {
    if(!this->Image->GetPixel(pixels[pixelId])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
      
    itk::VariableLengthVector<float> normalizedPixel;
    PixelType pixel = this->Image->GetPixel(pixels[pixelId]);
    normalizedPixel.SetSize(numberOfComponents);
    for(unsigned int component = 0; component < numberOfComponents; component++)
      {
      unsigned int channel = channelsToUse[component];
      normalizedPixel[component] = (pixel[channel] - minimumOfChannels[channel])/
                                   (maximumOfChannels[channel] - minimumOfChannels[channel]);
      if(this->Debug)
	{
	std::cout << "Pixel " << pixelId << " (" << pixels[pixelId] << ") channel " << channel << " has value " << pixel[channel] << " and normalized value " << normalizedPixel[component] << std::endl;
	debugNormalizedPixelValues.push_back(normalizedPixel[component]);
	}
      }
    
    sample->PushBack(normalizedPixel);
    }

  Helpers::WriteVectorToFile<float>(debugNormalizedPixelValues, "histogram.txt");
  
  histogramFilter->SetHistogramSize(histogramSize);
  histogramFilter->SetHistogramBinMinimum(binMinimum);
  histogramFilter->SetHistogramBinMaximum(binMaximum);
  histogramFilter->SetAutoMinimumMaximum(false);
  histogramFilter->SetInput(sample);
  histogramFilter->Modified();
  histogramFilter->Update();
  histogramFilter->Register();

  return histogramFilter->GetOutput();
  //histogramOutput = histogramFilter->GetOutput();
}

void ImageGraphCut::CreateHistograms()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms
  std::cout << "CreateHistograms()" << std::endl;
  
  std::vector<unsigned int> channelsToUse;
  if(this->IncludeColorInHistogram)
    {
    channelsToUse.push_back(0);
    channelsToUse.push_back(1);
    channelsToUse.push_back(2);
    }
  if(this->IncludeDepthInHistogram)
    {
    channelsToUse.push_back(3);
    }

  this->ForegroundHistogram = CreateHistogram(this->Sources, channelsToUse);
  this->BackgroundHistogram = CreateHistogram(this->Sinks, channelsToUse);
  
  //CreateHistogram(this->Sources, channelsToUse, this->ForegroundHistogram);
  //CreateHistogram(this->Sinks, channelsToUse, this->BackgroundHistogram);
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

  this->Sigma = ComputeAverageRandomDifferences(1000);
  
  if(this->Debug)
    {
    this->DebugGraphLines->Initialize();
  
    this->DebugGraphEdgeWeights->Initialize();
    this->DebugGraphEdgeWeights->SetNumberOfValues(1);
    }
  
  // We use a neighborhood iterator here even though we are looking only at a single pixel index in all images on each iteration because we use the neighborhood to determine edge validity.
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  NeighborhoodIteratorType iterator(Helpers::Get1x1Radius(), this->Image, this->Image->GetLargestPossibleRegion());
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
        float pixelDifference = this->DifferenceFunction->ComputeDifference(centerPixel, neighborPixel);

        // Compute the edge weight
        weight = ComputeNEdgeWeight(pixelDifference);

        }// end if current and neighbor are valid
      // Add the edge to the graph
      void* node1 = this->NodeImage->GetPixel(iterator.GetIndex());
      void* node2 = this->NodeImage->GetPixel(iterator.GetIndex(neighbors[i]));
      this->Graph->add_edge(node1, node2, weight, weight); // This is an undirected graph so we create a bidirectional edge with both weights set to 'weight'
      
      if(this->Debug)
	{
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0,this->DebugGraphPointIds->GetPixel(iterator.GetIndex()));
	line->GetPointIds()->SetId(1,this->DebugGraphPointIds->GetPixel(iterator.GetIndex(neighbors[i])));
	this->DebugGraphLines->InsertNextCell(line);
      
	this->DebugGraphEdgeWeights->InsertNextValue(weight);
	}

      //std::cout << "Set n-edge weight to " << weight << std::endl;
      } // end loop over neighbors
    } // end iteration over entire image

}

void ImageGraphCut::CreateTWeights()
{
  std::cout << "CreateTWeights()" << std::endl;
  ////////// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node) //////////

  // Compute the histograms of the selected foreground and background pixels
    
  CreateHistograms();
  
  std::vector<unsigned int> channelsToUse;
  if(this->IncludeColorInHistogram)
    {
    channelsToUse.push_back(0);
    channelsToUse.push_back(1);
    channelsToUse.push_back(2);
    }
  if(this->IncludeDepthInHistogram)
    {
    channelsToUse.push_back(3);
    }
      
  if(this->Debug)
    {
    std::cout << "Using channels ";
    for(unsigned int i = 0; i < channelsToUse.size(); ++i)
      {
      std::cout << channelsToUse[i] << " ";
      }
    std::cout << " to create T-Weights." << std::endl;
    unsigned int numberOfTuples = this->Image->GetLargestPossibleRegion().GetSize()[0] * this->Image->GetLargestPossibleRegion().GetSize()[1];
    this->DebugGraphSinkWeights->SetNumberOfTuples(numberOfTuples);
  
    this->DebugGraphSourceWeights->SetNumberOfTuples(numberOfTuples);
    
    this->DebugGraphSourceHistogram->SetNumberOfTuples(numberOfTuples);
    
    this->DebugGraphSinkHistogram->SetNumberOfTuples(numberOfTuples);

    }
  itk::ImageRegionIterator<ImageType> imageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NodeImageType> nodeIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  nodeIterator.GoToBegin();

  // Since the t-weight function takes the log of the histogram value,
  // we must handle bins with frequency = 0 specially (because log(0) = -inf)
  // For empty histogram bins we use tinyValue instead of 0.
  float tinyValue = 1e-10;
  
  // These are only for debuging/tracking
  std::vector<float> sinkTWeights;
  std::vector<float> sourceTWeights;
  std::vector<float> sourceHistogramValues;
  std::vector<float> sinkHistogramValues;

  std::vector<float> minimumOfChannels = Helpers::ComputeMinOfAllChannels(this->Image);
  std::vector<float> maximumOfChannels = Helpers::ComputeMaxOfAllChannels(this->Image);
  
  // Use the colors only for the t-weights
  unsigned int debugIteratorCounter = 0;
  while(!imageIterator.IsAtEnd())
    {
    PixelType pixel = imageIterator.Get();
    //float sinkHistogramValue = 0.0;
    //float sourceHistogramValue = 0.0;
    float sinkHistogramValue = tinyValue;
    float sourceHistogramValue = tinyValue;

    if(pixel[4]) // Pixel is valid
      {
      //std::cout << "Pixels have size: " << pixel.Size() << std::endl;

      HistogramType::MeasurementVectorType measurementVector(channelsToUse.size());
      for(unsigned int component = 0; component < channelsToUse.size(); component++)
        {
        unsigned int channel = channelsToUse[component];
        //measurementVector[component] = pixel[channel]; // Un-normalized

        measurementVector[component] = (pixel[channel] - minimumOfChannels[channel])/(maximumOfChannels[channel] - minimumOfChannels[channel]);
        }

      sinkHistogramValue = this->BackgroundHistogram->GetFrequency(this->BackgroundHistogram->GetIndex(measurementVector));
      sourceHistogramValue = this->ForegroundHistogram->GetFrequency(this->ForegroundHistogram->GetIndex(measurementVector));

      // Convert the histogram value/frequency to make it as if it came from a normalized histogram
      float normalizedSinkHistogramValue = sinkHistogramValue / static_cast<float>(this->BackgroundHistogram->GetTotalFrequency());
      float normalizedSourceHistogramValue = sourceHistogramValue / static_cast<float>(this->ForegroundHistogram->GetTotalFrequency());

      if(normalizedSinkHistogramValue <= 0)
        {
        normalizedSinkHistogramValue = tinyValue;
        }
      if(normalizedSourceHistogramValue <= 0)
        {
        normalizedSourceHistogramValue = tinyValue;
        }
	
//       std::cout << "Original value: " << pixel[3] << " normalized value: " << measurementVector[0] 
// 		<< " normalized source histogram count: " << normalizedSourceHistogramValue
// 		<< " normalized sink histogram count: " << normalizedSinkHistogramValue << std::endl;
// 		
      //std::cout << "Setting background weight to: " << -this->Lambda*log(sinkHistogramValue) << std::endl;
      //std::cout << "Setting foreground weight to: " << -this->Lambda*log(sourceHistogramValue) << std::endl;

      sinkHistogramValues.push_back(normalizedSinkHistogramValue);
      sourceHistogramValues.push_back(normalizedSourceHistogramValue);

      //float sinkWeight = -this->Lambda*log(normalizedSinkHistogramValue);
      // NOTE! The sink weight t-link is set as a function of the FOREGROUND probability.
      float sinkWeight = ComputeTEdgeWeight(Helpers::NegativeLog(normalizedSourceHistogramValue));
      sinkTWeights.push_back(sinkWeight);

      //float sourceWeight = -this->Lambda*log(normalizedSourceHistogramValue);
      // NOTE! The source weight t-link is set as a function of the BACKGROUND probability.
      float sourceWeight = ComputeTEdgeWeight(Helpers::NegativeLog(normalizedSinkHistogramValue));
      sourceTWeights.push_back(sourceWeight);

      // Add the edge to the graph and set its weight
      // See the table on p108 of "Interactive Graph Cuts for Optimal Boundary & Region Segmentation of Objects in N-D Images". 

      this->Graph->add_tweights(nodeIterator.Get(), sourceWeight, sinkWeight); // (node_id, source, sink)

      if(this->Debug)
        {
        this->DebugGraphSinkWeights->SetValue(debugIteratorCounter, sinkWeight);
        this->DebugGraphSourceWeights->SetValue(debugIteratorCounter, sourceWeight);

        this->DebugGraphSourceHistogram->SetValue(debugIteratorCounter, normalizedSourceHistogramValue);
        this->DebugGraphSinkHistogram->SetValue(debugIteratorCounter, normalizedSinkHistogramValue);
        }
      }
    else
      {
      this->Graph->add_tweights(nodeIterator.Get(), 0, 0);
      if(this->Debug)
        {
        this->DebugGraphSinkWeights->SetValue(debugIteratorCounter, 0);
        this->DebugGraphSourceWeights->SetValue(debugIteratorCounter, 0);
        this->DebugGraphSourceHistogram->SetValue(debugIteratorCounter, 0);
        this->DebugGraphSinkHistogram->SetValue(debugIteratorCounter, 0);
        }
      }
    debugIteratorCounter++;
    ++imageIterator;
    ++nodeIterator;
    }

  if(this->Debug)
    {
    std::cout << "Average sinkHistogramValue: " << Helpers::VectorAverage<float>(sinkHistogramValues) << std::endl;
    std::cout << "Average sourceHistogramValue: " << Helpers::VectorAverage<float>(sourceHistogramValues) << std::endl;
    
    std::cout << "Max sinkHistogramValue: " << *(std::max_element(sinkHistogramValues.begin(), sinkHistogramValues.end())) << std::endl;
    std::cout << "Max sourceHistogramValue: " << *(std::max_element(sourceHistogramValues.begin(), sourceHistogramValues.end())) << std::endl;
    
    std::cout << "Average sourceTWeights: " << Helpers::VectorAverage<float>(sourceTWeights) << std::endl;
    std::cout << "Average sinkTWeights: " << Helpers::VectorAverage<float>(sinkTWeights) << std::endl;
    
    std::cout << "Max sourceTWeights: " << *(std::max_element(sourceTWeights.begin(), sourceTWeights.end())) << std::endl;
    std::cout << "Max sinkTWeights: " << *(std::max_element(sinkTWeights.begin(), sinkTWeights.end())) << std::endl;
    }
}

void ImageGraphCut::SetHardSources(const std::vector<itk::Index<2> >& pixels)
{
  // Set very high source weights for the pixels which were selected as foreground by the user
  
  // If we are creating the debugging PolyData, we want to use the max of the "normal" t-weights instead of the infinity value so the range for visualization is reasonable
  float valuesRange[2];
  this->DebugGraphSourceWeights->GetValueRange(valuesRange);
  
  float highValue = std::numeric_limits<float>::max();
  //float highValue = 2.;
  // See the table on p108 of "Interactive Graph Cuts for Optimal Boundary & Region Segmentation of Objects in N-D Images". 
  // We want to set the source link high and the sink link to zero
  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    //std::cout << "Setting t-weight for node: " << this->NodeImage->GetPixel(pixels[i]) << " (pixel " << pixels[i] << ")" << std::endl;
    this->Graph->add_tweights(this->NodeImage->GetPixel(pixels[i]), highValue, 0); // (node_id, source, sink);
    }
}

void ImageGraphCut::SetHardSinks(const std::vector<itk::Index<2> >& pixels)
{
  // Set very high sink weights for the pixels which were selected as background by the user
  
  // See the table on p108 of "Interactive Graph Cuts for Optimal Boundary & Region Segmentation of Objects in N-D Images". 
  // We want to set the sink link high and the source link to zero. This means it is hard to cut the sink link, which is what we want.
  
  float highValue = std::numeric_limits<float>::max();
  //float highValue = 2.;
  
  // If we are creating the debugging PolyData, we want to use the max of the "normal" t-weights instead of the infinity value so the range for visualization is reasonable
  float valuesRange[2];
  this->DebugGraphSinkWeights->GetValueRange(valuesRange);
  
  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    this->Graph->add_tweights(this->NodeImage->GetPixel(pixels[i]), 0, highValue); // (node_id, source, sink);
    }
}

void ImageGraphCut::CreateGraph()
{
  if(this->Debug)
    {
    std::cout << "CreateGraph()" << std::endl;
    CreateDebugPolyData();
    }
    
  CreateGraphNodes();

  CreateNWeights();
  
  CreateTWeights();
  
  // Set very high source weights for the pixels which were selected as foreground by the user.
  SetHardSinks(this->Sinks);
  SetHardSources(this->Sources);

  if(this->Debug)
    {
    AssembleAndWriteDebugGraph();
    }
}

void ImageGraphCut::AssembleAndWriteDebugGraph()
{
  this->DebugGraphPolyData->SetLines(this->DebugGraphLines);
  this->DebugGraphPolyData->GetCellData()->SetScalars(this->DebugGraphEdgeWeights);
  this->DebugGraphPolyData->GetPointData()->AddArray(this->DebugGraphSinkWeights);
  this->DebugGraphPolyData->GetPointData()->AddArray(this->DebugGraphSourceWeights);
  this->DebugGraphPolyData->GetPointData()->AddArray(this->DebugGraphSourceHistogram);
  this->DebugGraphPolyData->GetPointData()->AddArray(this->DebugGraphSinkHistogram);
  
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("DebugGraph.vtp");
  writer->SetInputData(this->DebugGraphPolyData);
  writer->Write();
}

std::vector<itk::Index<2> > ImageGraphCut::GetSources()
{
  return this->Sources;
}

void ImageGraphCut::SetLambda(const float lambda)
{
  this->Lambda = lambda;
}

void ImageGraphCut::SetNumberOfHistogramBins(const int bins)
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

void ImageGraphCut::SetSources(vtkPolyData* const sources)
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

void ImageGraphCut::SetSinks(vtkPolyData* const sinks)
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

void ImageGraphCut::SetSources(const std::vector<itk::Index<2> >& sources)
{
  this->Sources = sources;
}

void ImageGraphCut::SetSinks(const std::vector<itk::Index<2> >& sinks)
{
  this->Sinks = sinks;
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


float ImageGraphCut::ComputeNEdgeWeight(const float difference)
{
  // This value should correspond to the variance (aka average) of the difference function you are using over the whole image.
  //float sigma = this->DifferenceFunction->AverageDifference;
  //float sigma = 1.0f;

  float sigma = this->Sigma;
  
  return exp(-pow(difference,2)/(2.0*sigma*sigma));
}

float ImageGraphCut::ComputeTEdgeWeight(const float value)
{
  return this->Lambda * value;
}

float ImageGraphCut::ComputeAverageRandomDifferences(const unsigned int numberOfDifferences)
{
  float sum = 0.0f;
  for(unsigned int i = 0; i < numberOfDifferences; ++i)
  {
    // Choose a random pixel
    itk::Index<2> pixel;
    pixel[0] = rand() % (this->Image->GetLargestPossibleRegion().GetSize()[0] - 2);
    pixel[1] = rand() % (this->Image->GetLargestPossibleRegion().GetSize()[1] - 2);

    itk::Index<2> pixelB = pixel;
    pixelB[0] += 1;

    if(!this->Image->GetLargestPossibleRegion().IsInside(pixel) || !this->Image->GetLargestPossibleRegion().IsInside(pixelB))
    {
      std::cout << "Pixel: " << pixel << " PixelB: " << pixelB << std::endl;
      std::cout << "Image: " << this->Image->GetLargestPossibleRegion() << std::endl;
      throw std::runtime_error("Something is wrong, pixels are not inside image!");
    }

    float difference = this->DifferenceFunction->ComputeDifference(this->Image->GetPixel(pixel), this->Image->GetPixel(pixelB));

    sum += difference;
  }

  return sum / static_cast<float>(numberOfDifferences);
}
