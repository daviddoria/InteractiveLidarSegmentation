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

#ifndef IMAGEGRAPHCUT_H
#define IMAGEGRAPHCUT_H

// VTK
#include <vtkSmartPointer.h>
class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;

// ITK
#include "itkImage.h"
#include "itkSampleToHistogramFilter.h"
#include "itkHistogram.h"
#include "itkListSample.h"

// STL
#include <vector>

// Custom
#include "Types.h"
#include "Difference.hpp"

// Kolmogorov's code
#include "graph.h"
typedef Graph GraphType;

// This is a special type to keep track of the graph node labels
typedef itk::Image<void*, 2> NodeImageType;

typedef itk::Statistics::Histogram< float,
        itk::Statistics::DenseFrequencyContainer2 > HistogramType;


class ImageGraphCut
{
public:
  ImageGraphCut();

  Difference* DifferenceFunction;
    
  float AverageDepthDifference;
  float MedianDepthDifference;
  float MinDepthDifference;
  float MaxDepthDifference;
  
  float AverageColorDifference;
  float MedianColorDifference;
  float MinColorDifference;
  float MaxColorDifference;
  
  // This function is purely for debugging. It lets you visualize the edge weights that will be used in the graph cut computation.
    
  // Several initializations are done here
  void SetImage(ImageType::Pointer image);
  ImageType::Pointer GetImage();
  
  // Create and cut the graph (The main driver function)
  void PerformSegmentation();

  // Get the masked output image
  ImageType::Pointer GetMaskedOutput();

  // Return a list of the selected (via scribbling) pixels
  std::vector<itk::Index<2> > GetSources();
  std::vector<itk::Index<2> > GetSinks();

  // Set the selected (via scribbling) pixels
  void SetSources(vtkPolyData* sources);
  void SetSinks(vtkPolyData* sinks);

  void SetSources(std::vector<itk::Index<2> > sources);
  void SetSinks(std::vector<itk::Index<2> > sinks);

  void SetHardSources(const std::vector<itk::Index<2> >& pixels);
  void SetHardSinks(const std::vector<itk::Index<2> >& pixels);
  
  // Get the output of the segmentation
  MaskImageType::Pointer GetSegmentMask();

  // Set the weight between the regional and boundary terms
  void SetLambda(float);

  // Set the number of bins per dimension of the foreground and background histograms
  void SetNumberOfHistogramBins(int);

  bool Debug;

  bool IncludeDepthInHistogram;
  bool IncludeColorInHistogram;
  
  bool IncludeDepthInDifference;
  bool IncludeColorInDifference;
  
  unsigned int NumberOfHistogramComponents;
  
  bool SecondStep;
  
  float BackgroundThreshold;

protected:

  float ComputeAverageRandomDifferences(const unsigned int numberOfDifferences);
  
  void CreateGraphNodes();
  
  // A Kolmogorov graph object
  GraphType* Graph;

  // The output segmentation
  MaskImageType::Pointer SegmentMask;

  // User specified foreground points
  std::vector<itk::Index<2> > Sources;

  // User specified background points
  std::vector<itk::Index<2> > Sinks;

  // The weighting between unary and binary terms
  float Lambda;

  // The number of bins per dimension of the foreground and background histograms
  int NumberOfHistogramBins;

  // An image which keeps tracks of the mapping between pixel index and graph node id
  NodeImageType::Pointer NodeImage;

  // Create the histograms from the users selections
  void CreateHistograms();
  const HistogramType* CreateHistogram(std::vector<itk::Index<2> > pixels, std::vector<unsigned int> channelsToUse);
  //void CreateHistogram(std::vector<itk::Index<2> > pixels, std::vector<unsigned int> channelsToUse, const HistogramType*);

  // Create a Kolmogorov graph structure from the image and selections
  void CreateGraphManually();
  void CreateGraph();
  void CreateNWeights();
  void CreateTWeights();
  
  // Perform the s-t min cut
  void CutGraph();

  // Several times throughout the algorithm we will need to traverse the image, looking exactly once at each edge. This iterator
  // creation is lengthy, so we do it once in this function and call it from everywhere we need it.
  void ConstructNeighborhoodIterator(NeighborhoodIteratorType* iterator, std::vector<NeighborhoodIteratorType::OffsetType>& neighbors);
  
  // Member variables

  // The histograms of the source and sink pixels
  const HistogramType* ForegroundHistogram;
  const HistogramType* BackgroundHistogram;

  // The image to be segmented
  ImageType::Pointer Image;

  // We need a 1x1 radius quite often, so we can create it in a single line with this function
  itk::Size<2> Get1x1Radius();
  
  // This function performs the negative exponential weighting
  float ComputeNEdgeWeight(const float difference);
  
  float ComputeTEdgeWeight(const float histogramValue);

  // Debugging variables/functions
  vtkSmartPointer<vtkPolyData> DebugGraphPolyData;
  vtkSmartPointer<vtkCellArray> DebugGraphLines;
  vtkSmartPointer<vtkFloatArray> DebugGraphEdgeWeights;
  vtkSmartPointer<vtkFloatArray> DebugGraphSourceWeights;
  vtkSmartPointer<vtkFloatArray> DebugGraphSinkWeights;
  vtkSmartPointer<vtkFloatArray> DebugGraphSourceHistogram;
  vtkSmartPointer<vtkFloatArray> DebugGraphSinkHistogram;
  itk::Image<unsigned int, 2>::Pointer DebugGraphPointIds;
  void CreateDebugPolyData();
  void AssembleAndWriteDebugGraph();
    
  std::vector<float> AllDepthDifferences;
  std::vector<float> AllColorDifferences;

  float Sigma;
};

#endif
