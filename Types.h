#ifndef TYPES_H
#define TYPES_H

#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImage.h"
//#include "itkShapedNeighborhoodIterator.h"
#include "itkVectorImage.h"

// All images are stored internally as float pixels
typedef itk::VectorImage<float,2> ImageType;
typedef itk::VariableLengthVector<float> PixelType;

typedef itk::Image<float, 2> FloatScalarImageType;

// For writing images, we need to first convert to actual unsigned char images
typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;
typedef UnsignedCharScalarImageType MaskImageType;

// For traversing image with an 8-neighborhood-visit-only-once idea
//typedef itk::ShapedNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
typedef itk::ConstShapedNeighborhoodIterator<ImageType> NeighborhoodIteratorType;

#endif