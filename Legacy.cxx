
void ImageGraphCut::CreateGraphManually()
{
  CreateGraphNodes();
  
  ////////// Create n-edges and set n-edge weights (links between image nodes) //////////
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
	//std::cout << "pixelDifference " << pixelDifference << std::endl;
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
      //std::cout << "Set n-edge weight to " << weight << std::endl;
      } // end loop over neighbors
    } // end iteration over entire image

  ////////// Add t-edges and set t-edge weights (links from image nodes to virtual background and virtual foreground node) //////////

  // Compute the histograms of the selected foreground and background pixels
  //
  //CreateFullHistogramSamples();
  if(this->IncludeDepthInHistogram)
    {
    CreateColorAndDepthHistogramSamples();
    }
  else
    {
    CreateColorHistogramSamples();
    }
    
  itk::ImageRegionIterator<ImageType> imageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NodeImageType> nodeIterator(this->NodeImage, this->NodeImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  nodeIterator.GoToBegin();

  // Since the t-weight function takes the log of the histogram value,
  // we must handle bins with frequency = 0 specially (because log(0) = -inf)
  // For empty histogram bins we use tinyValue instead of 0.
  float tinyValue = 1e-10;

  std::cout << "Setting T-Weights..." << std::endl;
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
      //unsigned int measurementVectorSize = std::min(pixel.Size(), 3u);
      unsigned int measurementVectorSize = 0;
      if(this->IncludeDepthInHistogram)
	{
	measurementVectorSize = 4; // This must be used if CreateFullHistogramSamples() is used
	}
      else
	{
	measurementVectorSize = 3; // This must be used if CreateColorHistogramSamples() is used
	}
      
      HistogramType::MeasurementVectorType measurementVector(measurementVectorSize);
      for(unsigned int i = 0; i < measurementVectorSize; i++)
	{
	//measurementVector[i] = pixel[i];
	measurementVector[i] = (pixel[i] - this->MinimumOfChannels[i])/(this->MaximumOfChannels[i] - this->MinimumOfChannels[i]);
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


void Difference::CreateNormalizedRGBImage(Vector3ImageType::Pointer image)
{
  image->SetRegions(this->Image->GetLargestPossibleRegion());
  image->Allocate();
  
  Vector3ImageType::PixelType zeroPixel;
  zeroPixel.Fill(0);
  
  image->FillBuffer(zeroPixel);
  
  itk::ImageRegionIterator<ImageType> fullImageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Vector3ImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  
  while(!fullImageIterator.IsAtEnd())
    {
    ImageType::PixelType fullPixel = fullImageIterator.Get();
  
    Vector3ImageType::PixelType rgbPixel;
    for(unsigned int i = 0; i < 3; ++i)
      {
      rgbPixel[i] = (fullPixel[i] - this->MinimumOfChannels[i]) / (this->MaximumOfChannels[i] - this->MinimumOfChannels[i]);
      }
  
    rgbImageIterator.Set(rgbPixel);
    ++fullImageIterator;
    ++rgbImageIterator;
    }
}


void Difference::CreateNormalizedDepthImage(FloatScalarImageType::Pointer image)
{
  image->SetRegions(this->Image->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);
  
  itk::ImageRegionIterator<ImageType> fullImageIterator(this->Image, this->Image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatScalarImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  
  while(!fullImageIterator.IsAtEnd())
    {
    ImageType::PixelType fullPixel = fullImageIterator.Get();
  
    float depthPixel = (fullPixel[3] - this->MinimumOfChannels[3]) / (this->MaximumOfChannels[3] - this->MinimumOfChannels[3]);
  
    rgbImageIterator.Set(depthPixel);
    ++fullImageIterator;
    ++rgbImageIterator;
    }
}



void ImageGraphCut::WriteEdgesManual(const std::string& fileName)
{
  std::cout << "WriteEdges()" << std::endl;
  
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
    p[2] = 0; // All points should have coordinate z=0 (lie on the XY plane)
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
  
  vtkSmartPointer<vtkFloatArray> graphEdgeWeights = vtkSmartPointer<vtkFloatArray>::New();
  graphEdgeWeights->SetNumberOfComponents(1);
  graphEdgeWeights->SetName("GraphEdgeWeights");
  
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
	
	graphEdgeWeights->InsertNextValue(ComputeEdgeWeight(depthDifference));
	
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
	graphEdgeWeights->InsertNextValue(0);
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
  polyData->GetCellData()->AddArray(graphEdgeWeights);
  
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(polyData);
  writer->Write();
}



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


double ImageGraphCut::ComputeNoise()
{
  // Compute an estimate of the "camera noise". This is used in the N-weight function.

  // This function should be idential to this->AverageColorDifference / NumberValidPixels
  
  // This function currently does not work !!!
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
      //DifferenceColorDataNormalized differenceColor(*(this->DifferenceFunction));
      /*
      DifferenceEuclidean differenceColor(*(this->DifferenceFunction));
      float colorDifference = differenceColor.Compute(centerPixel, neighborPixel);
      sigma += colorDifference;
      */
      numberOfEdges++;
      }

    }

  // Normalize
  sigma /= static_cast<double>(numberOfEdges);

  return sigma;
}

void ImageGraphCut::ComputeGlobalStatistics()
{
  ComputeMinAndMaxInAllChannels();
  
  std::cout << "Computed mins and maxs of all channels" << std::endl;

  
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


void ImageGraphCut::CreateColorHistogramSamples()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms

  // We want the histogram bins to take values from 0 to 255 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(3);
  HistogramType::MeasurementVectorType binMaximum(3);
  for(unsigned int i = 0; i < 3; i++)
    {
    binMinimum[i] = 0;
    binMaximum[i] = 255;
    }

  // Setup the histogram size
  //std::cout << "Image components per pixel: " << this->Image->GetNumberOfComponentsPerPixel() << std::endl;
  SampleToHistogramFilterType::HistogramSizeType histogramSize(3);
  histogramSize.Fill(this->NumberOfHistogramBins);

  // Create foreground samples and histogram
  this->ForegroundSample->Clear();
  this->ForegroundSample->SetMeasurementVectorSize(3);
  //std::cout << "Measurement vector size: " << this->ForegroundSample->GetMeasurementVectorSize() << std::endl;
  //std::cout << "Pixel size: " << this->Image->GetPixel(this->Sources[0]).GetNumberOfElements() << std::endl;
  
  for(unsigned int sourceCounter = 0; sourceCounter < this->Sources.size(); sourceCounter++)
    {
    if(!this->Image->GetPixel(this->Sources[sourceCounter])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
    itk::VariableLengthVector<float> pixel;
    pixel.SetSize(3);
    for(unsigned int component = 0; component < 3; component++)
      {
      pixel[component] = this->Image->GetPixel(this->Sources[sourceCounter])[component];
      }
    this->ForegroundSample->PushBack(pixel);
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
  this->BackgroundSample->SetMeasurementVectorSize(3);
  for(unsigned int sinkCounter = 0; sinkCounter < this->Sinks.size(); sinkCounter++)
    {
    if(!this->Image->GetPixel(this->Sinks[sinkCounter])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
    itk::VariableLengthVector<float> pixel;
    pixel.SetSize(3);
    for(unsigned int component = 0; component < 3; ++component)
      {
      pixel[component] = this->Image->GetPixel(this->Sinks[sinkCounter])[component];
      }
    this->BackgroundSample->PushBack(pixel);
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


void ImageGraphCut::CreateFullHistogramSamples()
{
  // This function creates ITK samples from the scribbled pixels and then computes the foreground and background histograms
  // This function doesn't make sense with RGBDV images because the validity channel should not be included in the histogram
  
  // We want the histogram bins to take values from 0 to 255 in all dimensions
  HistogramType::MeasurementVectorType binMinimum(this->Image->GetNumberOfComponentsPerPixel());
  HistogramType::MeasurementVectorType binMaximum(this->Image->GetNumberOfComponentsPerPixel());
  for(unsigned int i = 0; i < this->Image->GetNumberOfComponentsPerPixel(); i++)
    {
    binMinimum[i] = 0; // These should instead be normalized [0,1]!!!
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
    if(!this->Image->GetPixel(this->Sources[i])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
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
    if(!this->Image->GetPixel(this->Sinks[i])[4]) // Don't include invalid pixels in the histogram
      {
      continue;
      }
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


/*
 // Display segmented image with black background pixels
void Form::StopProgressSlot()
{
  // When the ProgressThread emits the StopProgressSignal, we need to display the result of the segmentation

  // Convert the masked image into a VTK image for display
  vtkSmartPointer<vtkImageData> VTKSegmentMask =
    vtkSmartPointer<vtkImageData>::New();
  if(this->GraphCut->GetPixelDimensionality() == 1)
    {
    ITKImagetoVTKImage<GrayscaleImageType>(static_cast<ImageGraphCut<GrayscaleImageType>* >(this->GraphCut)->GetMaskedOutput(), VTKSegmentMask);
    }
  else if(this->GraphCut->GetPixelDimensionality() == 3)
    {
    ITKImagetoVTKImage<ColorImageType>(static_cast<ImageGraphCut<ColorImageType>* >(this->GraphCut)->GetMaskedOutput(), VTKSegmentMask);
    }
  else if(this->GraphCut->GetPixelDimensionality() == 5)
    {
    ITKImagetoVTKImage<RGBDIImageType>(static_cast<ImageGraphCut<RGBDIImageType>* >(this->GraphCut)->GetMaskedOutput(), VTKSegmentMask);
    }
  else
    {
    std::cerr << "This type of image (" << this->GraphCut->GetPixelDimensionality() << ") cannot be displayed!" << std::endl;
    exit(-1);
    }

  // Remove the old output, set the new output and refresh everything
  this->ResultActor = vtkSmartPointer<vtkImageActor>::New();
  this->ResultActor->SetInput(VTKSegmentMask);
  this->RightRenderer->RemoveAllViewProps();
  this->RightRenderer->AddActor(ResultActor);
  this->RightRenderer->ResetCamera();
  this->Refresh();

  this->progressBar->hide();
}
*/


void MainWindow::on_actionSaveSegmentation_triggered()
{
  // Ask the user for a filename to save the segment mask image to

  QString fileName = QFileDialog::getSaveFileName(this,
    "Save Segment Mask Image", ".", "PNG Files (*.png)");
/*
  // Convert the image from a 1D vector image to an unsigned char image
  typedef itk::CastImageFilter< GrayscaleImageType, itk::Image<itk::CovariantVector<unsigned char, 1>, 2 > > CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(this->GraphCut->GetSegmentMask());

  typedef itk::NthElementImageAdaptor< itk::Image<itk:: CovariantVector<unsigned char, 1>, 2 >,
    unsigned char> ImageAdaptorType;

  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SelectNthElement(0);
  adaptor->SetImage(castFilter->GetOutput());
*/
/*
  // Write the file (object is white)
  //typedef  itk::ImageFileWriter< ImageAdaptorType > WriterType;
  typedef  itk::ImageFileWriter<Mask> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  //writer->SetInput(adaptor);
  writer->SetInput(this->GraphCut.GetSegmentMask());
  writer->Update();
  */

  // Write the inverted file (object is black)
  Mask::Pointer inverted = Mask::New();
  Helpers::InvertBinaryImage(this->GraphCut.GetSegmentMask(), inverted);

  //typedef  itk::ImageFileWriter< ImageAdaptorType > WriterType;
  typedef  itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  writer->SetInput(inverted);
  writer->Update();

}


void LidarSegmentationWidget::GenerateNeighborSinks()
{
  MaskImageType::Pointer sourcesImage = MaskImageType::New();
  sourcesImage->SetRegions(this->ImageRegion);
  Helpers::IndicesToBinaryImage(this->Sources, sourcesImage);

  // Dilate the mask
  std::cout << "Dilating mask..." << std::endl;
  typedef itk::BinaryBallStructuringElement<MaskImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(1);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <MaskImageType, MaskImageType, StructuringElementType> BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(sourcesImage);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  // Binary XOR the images to get the difference image
  std::cout << "XORing masks..." << std::endl;
  typedef itk::XorImageFilter <MaskImageType> XorImageFilterType;

  XorImageFilterType::Pointer xorFilter = XorImageFilterType::New();
  xorFilter->SetInput1(dilateFilter->GetOutput());
  xorFilter->SetInput2(sourcesImage);
  xorFilter->Update();

  Helpers::WriteImage<MaskImageType>(xorFilter->GetOutput(), "boundaryOfSegmentation.png");

  /*
  // Dilate the mask
  std::cout << "Dilating mask..." << std::endl;
  typedef itk::BinaryBallStructuringElement<MaskImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElementBig;
  structuringElementBig.SetRadius(7);
  structuringElementBig.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <MaskImageType, MaskImageType, StructuringElementType> BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilterBig = BinaryDilateImageFilterType::New();
  dilateFilterBig->SetInput(sourcesImage);
  dilateFilterBig->SetKernel(structuringElementBig);
  dilateFilterBig->Update();

  StructuringElementType structuringElementSmall;
  structuringElementSmall.SetRadius(6);
  structuringElementSmall.CreateStructuringElement();

  BinaryDilateImageFilterType::Pointer dilateFilterSmall = BinaryDilateImageFilterType::New();
  dilateFilterSmall->SetInput(sourcesImage);
  dilateFilterSmall->SetKernel(structuringElementSmall);
  dilateFilterSmall->Update();

  //Helpers::WriteImage<MaskImageType>(dilateFilter->GetOutput(), "dilated.png");

  // Binary XOR the images to get the difference image
  std::cout << "XORing masks..." << std::endl;
  typedef itk::XorImageFilter <MaskImageType> XorImageFilterType;

  XorImageFilterType::Pointer xorFilter = XorImageFilterType::New();
  xorFilter->SetInput1(dilateFilterBig->GetOutput());
  xorFilter->SetInput2(dilateFilterSmall->GetOutput());
  xorFilter->Update();

  Helpers::WriteImage<MaskImageType>(xorFilter->GetOutput(), "boundaryOfSegmentation.png");
  */

  // Iterate over the border pixels. If the closest pixel in the original segmentation has
  // a depth greater than a threshold, mark it as a new sink. Else, do not.
  std::cout << "Determining which boundary pixels should be declared background..." << std::endl;
  //std::cout << "There should be " << Helpers::CountNonZeroPixels(xorFilter->GetOutput())
  //          << " considered." << std::endl;

  /*
  std::vector<itk::Index<2> > newSinks;
  itk::ImageRegionIterator<MaskImageType> imageIterator(xorFilter->GetOutput(),
                                     xorFilter->GetOutput()->GetLargestPossibleRegion());

  unsigned int consideredCounter = 0;
  unsigned int backgroundCounter = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // If the current pixel is in question
      {
      consideredCounter++;
      //std::cout << "Considering pixel " << consideredCounter << " (index "
                  << imageIterator.GetIndex() << ")" << std::endl;
      ImageType::PixelType currentPixel = this->Image->GetPixel(imageIterator.GetIndex());
      itk::Index<2> closestPixelIndex = Helpers::FindClosestNonZeroPixel(sourcesImage, imageIterator.GetIndex());
      //std::cout << "Closest pixel is " << closestPixelIndex << std::endl;
      ImageType::PixelType closestPixel = this->Image->GetPixel(closestPixelIndex);
      //std::cout << "Current pixel depth value is " << currentPixel[3] << std::endl;
      //std::cout << "Closest pixel depth value is " << closestPixel[3] << std::endl;
      float difference = fabs(currentPixel[3]-closestPixel[3]);
      if(difference > this->txtBackgroundThreshold->text().toFloat())
        {
        //std::cout << "Difference was " << difference << " so this is a sink pixel." << std::endl;
        newSinks.push_back(imageIterator.GetIndex());
        backgroundCounter++;
        }
      else
        {
        //std::cout << "Difference was " << difference << " so this is NOT a sink pixel." << std::endl;
        }
      }

    ++imageIterator;
    }
  */
  std::vector<itk::Index<2> > newSinks;
  itk::ImageRegionIterator<MaskImageType> imageIterator(xorFilter->GetOutput(),
                                                        xorFilter->GetOutput()->GetLargestPossibleRegion());

  typedef itk::VectorIndexSelectionCastImageFilter<ImageType, FloatScalarImageType> IndexSelectionType;
  IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(3);
  indexSelectionFilter->SetInput(this->Image);
  indexSelectionFilter->Update();

  FloatScalarImageType::Pointer depthImage = indexSelectionFilter->GetOutput();

  std::vector<itk::Index<2> > consideredPixels;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // If the current pixel is in question
      {
      consideredPixels.push_back(imageIterator.GetIndex());
      //std::cout << "Considering pixel " << consideredCounter << " (index "
      //          << imageIterator.GetIndex() << ")" << std::endl;
      ImageType::PixelType currentPixel = this->Image->GetPixel(imageIterator.GetIndex());

      typedef itk::RegionOfInterestImageFilter< FloatScalarImageType, FloatScalarImageType >
              RegionOfInterestImageFilterType;
      RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
      //unsigned int radius = 2;
      unsigned int radius = this->txtBackgroundCheckRadius->text().toUInt();
      ImageType::IndexType start;
      start[0] = imageIterator.GetIndex()[0] - radius;
      start[1] = imageIterator.GetIndex()[1] - radius;

      ImageType::SizeType size;
      size.Fill(radius*2 + 1);

      ImageType::RegionType desiredRegion(start, size);
      //std::cout << "Checking pixel: " << imageIterator.GetIndex()
      //          << " with region: " << desiredRegion << std::endl;

      regionOfInterestImageFilter->SetRegionOfInterest(desiredRegion);
      regionOfInterestImageFilter->SetInput(depthImage);
      regionOfInterestImageFilter->Update();

      typedef itk::StatisticsImageFilter<FloatScalarImageType> StatisticsImageFilterType;
      StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
      statisticsImageFilter->SetInput(regionOfInterestImageFilter->GetOutput());
      statisticsImageFilter->Update();

      float difference = statisticsImageFilter->GetMaximum() - statisticsImageFilter->GetMinimum();
      if(this->Debug)
        {
        std::cout << "Max: " << statisticsImageFilter->GetMaximum()
                  << " min: " <<  statisticsImageFilter->GetMinimum() << std::endl;
        }

      if(difference > this->txtBackgroundThreshold->text().toFloat())
        {
        //std::cout << "Difference was " << difference << " so this is a sink pixel." << std::endl;
        newSinks.push_back(imageIterator.GetIndex());
        }
      else
        {
        //std::cout << "Difference was " << difference << " so this is NOT a sink pixel." << std::endl;
        }
      }

    ++imageIterator;
    }

  unsigned char blue[3] = {0, 0, 255};

  Helpers::SetPixels(this->SourceSinkImageData, consideredPixels, blue);
  this->SourceSinkImageData->Modified();
  this->Refresh();

  // Save the new sink pixels for inspection
  UnsignedCharScalarImageType::Pointer newSinksImage = UnsignedCharScalarImageType::New();
  newSinksImage->SetRegions(this->Image->GetLargestPossibleRegion());
  newSinksImage->Allocate();

  Helpers::IndicesToBinaryImage(newSinks, newSinksImage);
  Helpers::WriteImage<MaskImageType>(newSinksImage, "newSinks.png");

  //std::cout << "Out of " << consideredCounter << " pixels considered, "
  //          << backgroundCounter << " were declared background." << std::endl;
  // Set the new sinks
  std::cout << "Setting " << newSinks.size() << " new sinks." << std::endl;

  // Modify the list of sinks so it can be retrieved by the MainWindow after the segmentation is finished
  this->Sinks.insert(this->Sinks.end(), newSinks.begin(), newSinks.end());

  UpdateSelections();
}
