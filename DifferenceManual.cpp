

  // If the image is more than 3 channels, compute the weighted difference,
  // weighting the first 3 channels by RGB_Weight/3 and the remaining channels by (1-RGB_Weight)/(NumChannels-3)
  //float PixelDifference(PixelType, PixelType);

  // Compute the difference of the first 3 channels (assumed to be R, G, and B)
  //float PixelColorDifference(PixelType, PixelType);
  //float HSVColorDifference(PixelType, PixelType);
  //float CIELABColorDifference(PixelType, PixelType);

  // Compute the difference of the first 3 channels (assumed to be R, G, and B)
  //float DataNormalizedPixelColorDifference(PixelType, PixelType);
  //float AbsoluteNormalizedPixelColorDifference(PixelType, PixelType);
  
  // Compute the difference of the 4th channel (channel 3) (assumed to be Depth)
  //float PixelDepthDifference(PixelType, PixelType);
  
  // Compute the difference of the 4th channel (channel 3) (assumed to be Depth)
  //float AbsoluteNormalizedPixelDepthDifference(PixelType, PixelType);
  //float DataNormalizedPixelDepthDifference(PixelType, PixelType);

  // Compute the normalized color difference and normalized depth difference and return the largest.
  //float PixelDifferenceMaxOfColorOrDepth(PixelType, PixelType);
  
  // Compute the 
  //float PixelDifferenceDepthWeightedByColor(PixelType, PixelType);
  
  // Compute means and medians of color and depth differences
  void ComputeGlobalStatistics();
  
  // Compute the average color difference between all pairs of neighboring pixels
  float ComputeAverageColorDifference();

  // Compute the average depth difference between all pairs of neighboring pixels
  float ComputeAverageDepthDifference();

  // Compute the median color difference between all pairs of neighboring pixels
  float ComputeMedianColorDifference();

  // Compute the median depth difference between all pairs of neighboring pixels
  float ComputeMedianDepthDifference();
  
  // Compute all depth differences
  void ComputeAllDepthDifferences();
  
  // Compute all color differences
  void ComputeAllColorDifferences();

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
