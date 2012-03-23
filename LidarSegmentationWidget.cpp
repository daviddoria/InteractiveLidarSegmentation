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

#include "LidarSegmentationWidget.h"

// Custom
#include "Difference.hpp"
#include "Helpers.h"
#include "InteractorStyleImageNoLevel.h"
#include "InteractorStyleScribble.h"

// ITK
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLineIterator.h"
#include "itkNthElementImageAdaptor.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkXorImageFilter.h"

// VTK
#include <vtkCamera.h>
#include <vtkImageData.h>
#include <vtkImageProperty.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>

// Qt
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>
#include <QTimer>
#include <QtConcurrentRun>

// STL
#include <iostream>

LidarSegmentationWidget::LidarSegmentationWidget(QWidget *parent)
{
  SharedConstructor();
}

LidarSegmentationWidget::LidarSegmentationWidget(const std::string& fileName)
{
  SharedConstructor();
  OpenFile(fileName);
}

void LidarSegmentationWidget::SharedConstructor()
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);

  this->ProgressDialog = new QProgressDialog();
  this->ProgressDialog->setMinimum(0);
  this->ProgressDialog->setMaximum(0);
  this->ProgressDialog->setWindowModality(Qt::WindowModal);
  connect(&this->FutureWatcher, SIGNAL(finished()), this, SLOT(slot_SegmentationComplete()));
  connect(&this->FutureWatcher, SIGNAL(finished()), this->ProgressDialog , SLOT(cancel()));
  
  // Global settings
  this->Flipped = false;
  this->Debug = true;
  
  // Qt connections
  // connect( this->sldHistogramBins, SIGNAL( valueChanged(int) ), this, SLOT(sldHistogramBins_valueChanged()));
  connect( this->sldLambda, SIGNAL( valueChanged(int) ), this, SLOT(UpdateLambda()));
  connect( this->txtLambdaMax, SIGNAL( textEdited(QString) ), this, SLOT(UpdateLambda()));

  // Setup backgrounds
  this->BackgroundColor[0] = 0;
  this->BackgroundColor[1] = 0;
  this->BackgroundColor[2] = .5;

  // Left pane
  this->OriginalImageData = vtkSmartPointer<vtkImageData>::New();
  this->OriginalImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->OriginalImageSliceMapper->SetInputData(this->OriginalImageData);
  this->OriginalImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  // Make the pixels sharp instead of blurry when zoomed
  this->OriginalImageSlice->GetProperty()->SetInterpolationTypeToNearest();
  this->OriginalImageSlice->SetMapper(this->OriginalImageSliceMapper);

  this->LeftRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->LeftRenderer->GradientBackgroundOn();
  this->LeftRenderer->SetBackground(this->BackgroundColor);
  this->LeftRenderer->SetBackground2(1,1,1);
  this->qvtkWidgetLeft->GetRenderWindow()->AddRenderer(this->LeftRenderer);

  this->LeftRenderer->AddViewProp(this->OriginalImageSlice);

  this->LeftInteractorStyle = vtkSmartPointer<InteractorStyleScribble>::New();
  this->LeftInteractorStyle->AddObserver(this->LeftInteractorStyle->ScribbleEvent,
                                         this, &LidarSegmentationWidget::ScribbleEventHandler);
  this->LeftInteractorStyle->SetCurrentRenderer(this->LeftRenderer);
  this->qvtkWidgetLeft->GetInteractor()->SetInteractorStyle(this->LeftInteractorStyle);

  // Right pane
  this->ResultImageData = vtkSmartPointer<vtkImageData>::New();
  this->ResultImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->ResultImageSliceMapper->SetInputData(this->ResultImageData);
  this->ResultImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  this->ResultImageSlice->SetMapper(this->ResultImageSliceMapper);
  
  this->RightRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->RightRenderer->GradientBackgroundOn();
  this->RightRenderer->SetBackground(this->BackgroundColor);
  this->RightRenderer->SetBackground2(1,1,1);
  this->qvtkWidgetRight->GetRenderWindow()->AddRenderer(this->RightRenderer);
  
  this->RightRenderer->AddViewProp(this->ResultImageSlice);

  this->RightInteractorStyle = vtkSmartPointer<InteractorStyleImageNoLevel>::New();
  this->RightInteractorStyle->SetCurrentRenderer(this->RightRenderer);
  this->qvtkWidgetRight->GetInteractor()->SetInteractorStyle(this->RightInteractorStyle);

  // Both panes
  this->SourceSinkImageData = vtkSmartPointer<vtkImageData>::New();
  
  this->LeftSourceSinkImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->LeftSourceSinkImageSliceMapper->SetInputData(this->SourceSinkImageData);
  
  this->LeftSourceSinkImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  // Make the pixels sharp instead of blurry when zoomed
  this->LeftSourceSinkImageSlice->GetProperty()->SetInterpolationTypeToNearest(); 
  this->LeftSourceSinkImageSlice->SetMapper(this->LeftSourceSinkImageSliceMapper);
  
  this->LeftRenderer->AddViewProp(this->LeftSourceSinkImageSlice);
  
  this->RightSourceSinkImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->RightSourceSinkImageSliceMapper->SetInputData(this->SourceSinkImageData);
  
  this->RightSourceSinkImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  // Make the pixels sharp instead of blurry when zoomed
  this->RightSourceSinkImageSlice->GetProperty()->SetInterpolationTypeToNearest(); 
  this->RightSourceSinkImageSlice->SetMapper(this->RightSourceSinkImageSliceMapper);

  this->RightRenderer->AddViewProp(this->RightSourceSinkImageSlice);
  
  // Default GUI settings
  this->radForeground->setChecked(true);

  // Setup toolbar
  // Open file buttons
  QIcon openIcon = QIcon::fromTheme("document-open");
  actionOpenImage->setIcon(openIcon);
  this->toolBar->addAction(actionOpenImage);

  UpdateLambda();
  
  // Save buttons
  QIcon saveIcon = QIcon::fromTheme("document-save");
  actionSaveSegmentation->setIcon(saveIcon);
  this->toolBar->addAction(actionSaveSegmentation);

  // Initialize GUI components
  //sldHistogramBins_valueChanged();
  //QTimer::singleShot(0, this, SLOT(sldHistogramBins_valueChanged()));
  //QTimer::singleShot(0, this, SLOT(sldHistogramBins_sliderMoved()));
  this->lblHistogramBins->setNum(this->sldHistogramBins->value());
  
  this->Image = ImageType::New();
}

void LidarSegmentationWidget::on_actionExit_triggered()
{
  exit(0);
}

void LidarSegmentationWidget::SetCameraPosition1()
{
  double leftToRight[3] = {-1,0,0};
  double bottomToTop[3] = {0,1,0};
  this->LeftInteractorStyle->SetImageOrientation(leftToRight, bottomToTop);
  this->RightInteractorStyle->SetImageOrientation(leftToRight, bottomToTop); 
  
  this->LeftRenderer->ResetCamera();
  this->LeftRenderer->ResetCameraClippingRange();
  
  this->RightRenderer->ResetCamera();
  this->RightRenderer->ResetCameraClippingRange();
}

void LidarSegmentationWidget::SetCameraPosition2()
{
  double leftToRight[3] = {-1,0,0};
  double bottomToTop[3] = {0,-1,0};
  this->LeftInteractorStyle->SetImageOrientation(leftToRight, bottomToTop);
  this->RightInteractorStyle->SetImageOrientation(leftToRight, bottomToTop); 
  
  this->LeftRenderer->ResetCamera();
  this->LeftRenderer->ResetCameraClippingRange();
  
  this->RightRenderer->ResetCamera();
  this->RightRenderer->ResetCameraClippingRange();
}

void LidarSegmentationWidget::on_actionFlipImage_triggered()
{
  if(this->Flipped)
    {
    SetCameraPosition1();
    }
  else
    {
    SetCameraPosition2();
    }
    
  this->Flipped = !this->Flipped;
  
  this->Refresh();
}

void LidarSegmentationWidget::on_actionSaveSegmentation_triggered()
{
  // Ask the user for a filename to save the segment mask image to

  QString fileName = QFileDialog::getSaveFileName(this,
    "Save Segment Mask Image", ".", "PNG Files (*.png)");

  if(fileName.isEmpty())
  {
    return;
  }
  // Write the file (object is white)
  typedef  itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  writer->SetInput(this->GraphCut.GetSegmentMask());
  writer->Update();
  
  /*
  // Write the inverted file (object is black)
  MaskImageType::Pointer inverted = MaskImageType::New();
  Helpers::InvertBinaryImage(this->GraphCut.GetSegmentMask(), inverted);
  
  //typedef  itk::ImageFileWriter< ImageAdaptorType > WriterType;
  typedef  itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  writer->SetInput(inverted);
  writer->Update();
  */
}

void LidarSegmentationWidget::on_actionOpenImage_triggered()
{
  //std::cout << "actionOpenImage_triggered()" << std::endl;
  //std::cout << "Enter OpenFile()" << std::endl;

  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this,
     "Open Image", ".", "RGBD Files (*.mha)");

  if(filename.isEmpty())
    {
    return;
    }
    
  OpenFile(filename.toStdString());
}

void LidarSegmentationWidget::DisplaySegmentationResult()
{
  // Convert the segmentation mask to a binary VTK image
  vtkSmartPointer<vtkImageData> VTKSegmentMask = vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKScalarImageToVTKImage(this->GraphCut.GetSegmentMask(), VTKSegmentMask);

  // Convert the image into a VTK image for display
  vtkSmartPointer<vtkImageData> VTKImage = vtkSmartPointer<vtkImageData>::New();
  //Helpers::ITKImageToVTKImage(this->GraphCut.GetMaskedOutput(), VTKImage);
  
  ImageType::Pointer maskedImage = ImageType::New();
  Helpers::MaskImage(this->Image.GetPointer(), this->GraphCut.GetSegmentMask(), maskedImage.GetPointer());
  Helpers::ITKImageToVTKImage(maskedImage.GetPointer(), VTKImage);
  
  // Mask the VTK image with the segmentation result/mask
  vtkSmartPointer<vtkImageData> VTKMaskedImage = vtkSmartPointer<vtkImageData>::New();
  Helpers::MaskImage(VTKImage, VTKSegmentMask, VTKMaskedImage);

  // Set the new output and refresh everything
  this->ResultImageSliceMapper->SetInputData(VTKMaskedImage);
  
  // Sometimes (in the middle of a two-step segmentation) sources/sinks are modified by the GraphCut object
  this->Sources = this->GraphCut.GetSources();
  this->Sinks = this->GraphCut.GetSinks();
  UpdateSelections();
  
  // We currently remove everything and re-add everything - see note about the "should be unnecessary" lines below.
  //this->RightRenderer->RemoveAllViewProps();
  //this->RightRenderer->AddViewProp(this->RightSourceSinkImageSlice);
  
  // These should be unnecessary because they are connected in the constructor... but the output is not displayed without them
  //this->ResultImageSlice->SetMapper(this->ResultImageSliceMapper);
  //this->RightRenderer->AddViewProp(this->ResultImageSlice);

  
  //this->RightRenderer->ResetCamera();
  this->Refresh();
}

void LidarSegmentationWidget::ScribbleEventHandler(vtkObject* caller, long unsigned int eventId, void* callData)
{
  //std::cout << "Handled scribble event." << std::endl;

  // Only mark the pixels on the path as foreground pixels
  //std::vector<itk::Index<2> > selection = this->LeftInteractorStyle->GetSelection();

  // Dilate the path and mark the path and the dilated pixel of the path as foreground
  unsigned int dilateRadius = 2;
  std::vector<itk::Index<2> > thinSelection = this->LeftInteractorStyle->GetSelection();
  std::vector<itk::Index<2> > selection = Helpers::DilatePixelList(thinSelection,
                                                                   this->Image->GetLargestPossibleRegion(),
                                                                   dilateRadius);
  
  if(this->radForeground->isChecked())
    {
    this->Sources.insert(this->Sources.end(), selection.begin(), selection.end());
    }
  else if(this->radBackground->isChecked())
    {
    this->Sinks.insert(this->Sinks.end(), selection.begin(), selection.end());
    }
  else
    {
    throw std::runtime_error("Something is wrong - either Foreground or Background selection mode must be selected.");
    }
    
  UpdateSelections();
}

void LidarSegmentationWidget::UpdateSelections()
{
  // First, clear the image
  Helpers::CreateTransparentImage(this->SourceSinkImageData);
  
  unsigned char green[3] = {0, 255, 0};
  unsigned char red[3] = {255, 0, 0};
  
  Helpers::SetPixels(this->SourceSinkImageData, this->Sources, green);
  Helpers::SetPixels(this->SourceSinkImageData, this->Sinks, red);

  this->SourceSinkImageData->Modified();
  
  std::cout << this->Sources.size() << " sources." << std::endl;
  std::cout << this->Sinks.size() << " sinks." << std::endl;
  
  //this->LeftSourceSinkImageSliceMapper->Modified();
  //this->RightSourceSinkImageSliceMapper->Modified();
  
  this->Refresh();
}

void LidarSegmentationWidget::slot_SegmentationComplete()
{
  // Display the result of the segmentation
  DisplaySegmentationResult();
}

float LidarSegmentationWidget::ComputeLambda()
{
  // Compute lambda by multiplying the percentage set by the slider by the MaxLambda set in the text box

  double lambdaMax = this->txtLambdaMax->text().toDouble();
  double lambdaPercent = this->sldLambda->value()/100.;
  double lambda = lambdaPercent * lambdaMax;

  return lambda;
}

void LidarSegmentationWidget::UpdateLambda()
{
  // Compute lambda and then set the label to this value so the user can see the current setting
  double lambda = ComputeLambda();
  this->lblLambda->setText(QString::number(lambda));
}

void LidarSegmentationWidget::on_sldHistogramBins_valueChanged()
{
  this->GraphCut.SetNumberOfHistogramBins(sldHistogramBins->value());
  //this->lblHistogramBins->setText(QString::number(sldHistogramBins->value())); // This is taken care of by a signal/slot pair setup in QtDesigner
}

void LidarSegmentationWidget::on_radForeground_clicked()
{
  this->LeftInteractorStyle->SetColorToGreen();
}

void LidarSegmentationWidget::on_radBackground_clicked()
{
  this->LeftInteractorStyle->SetColorToRed();
}

void LidarSegmentationWidget::on_btnClearSelections_clicked()
{
  //this->LeftInteractorStyle->ClearSelections();
  this->Sources.clear();
  this->Sinks.clear();
  UpdateSelections();
}

void LidarSegmentationWidget::on_btnClearForeground_clicked()
{
  //this->LeftInteractorStyle->ClearForegroundSelections();
  this->Sources.clear();
  UpdateSelections();
}

void LidarSegmentationWidget::on_btnClearBackground_clicked()
{
  //this->LeftInteractorStyle->ClearBackgroundSelections();
  this->Sinks.clear();
  UpdateSelections();
}

void LidarSegmentationWidget::on_action_Selections_SaveSelectionsAsImage_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this,
     "Save Image", ".", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  RGBImageType::Pointer selectionsImage = RGBImageType::New();
  
  selectionsImage->SetRegions(this->ImageRegion);
  selectionsImage->Allocate();
  
  RGBPixelType whitePixel;
  whitePixel.SetRed(255);
  whitePixel.SetGreen(255);
  whitePixel.SetBlue(255);
  
  selectionsImage->FillBuffer(whitePixel);
  
  RGBPixelType greenPixel;
  greenPixel.SetRed(0);
  greenPixel.SetGreen(255);
  greenPixel.SetBlue(0);
  Helpers::SetPixels<RGBImageType>(selectionsImage, this->Sources, greenPixel);
  
  RGBPixelType redPixel;
  redPixel.SetRed(255);
  redPixel.SetGreen(0);
  redPixel.SetBlue(0);
  Helpers::SetPixels<RGBImageType>(selectionsImage, this->Sinks, redPixel);

  typedef  itk::ImageFileWriter< RGBImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.toStdString());
  writer->SetInput(selectionsImage);
  writer->Update();
}


void LidarSegmentationWidget::on_action_Selections_SaveForegroundSelectionsAsImage_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this,
     "Save Image", "foreground.png", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  RGBImageType::Pointer selectionsImage = RGBImageType::New();

  selectionsImage->SetRegions(this->ImageRegion);
  selectionsImage->Allocate();

  RGBPixelType blackPixel;
  blackPixel.SetRed(0);
  blackPixel.SetGreen(0);
  blackPixel.SetBlue(0);
  
  RGBPixelType whitePixel;
  whitePixel.SetRed(255);
  whitePixel.SetGreen(255);
  whitePixel.SetBlue(255);

  Helpers::SetPixelsInRegionToValue(selectionsImage.GetPointer(), selectionsImage->GetLargestPossibleRegion(),
                                    blackPixel);

  Helpers::SetPixels(selectionsImage.GetPointer(), this->Sources, whitePixel);

  typedef  itk::ImageFileWriter< RGBImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.toStdString());
  writer->SetInput(selectionsImage);
  writer->Update();
}


void LidarSegmentationWidget::on_action_Selections_SaveBackgroundSelectionsAsImage_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this,
     "Save Image", "background.png", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  RGBImageType::Pointer selectionsImage = RGBImageType::New();

  selectionsImage->SetRegions(this->ImageRegion);
  selectionsImage->Allocate();

  RGBPixelType blackPixel;
  blackPixel.SetRed(0);
  blackPixel.SetGreen(0);
  blackPixel.SetBlue(0);

  RGBPixelType whitePixel;
  whitePixel.SetRed(255);
  whitePixel.SetGreen(255);
  whitePixel.SetBlue(255);

  Helpers::SetPixelsInRegionToValue(selectionsImage.GetPointer(), selectionsImage->GetLargestPossibleRegion(),
                                    blackPixel);

  Helpers::SetPixels(selectionsImage.GetPointer(), this->Sinks, whitePixel);

  typedef  itk::ImageFileWriter< RGBImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.toStdString());
  writer->SetInput(selectionsImage);
  writer->Update();
}

void LidarSegmentationWidget::on_action_Selections_SaveSelectionsAsText_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this,
     "Save Selections", ".", "TXT Files (*.txt)");

  if(filename.isEmpty())
    {
    return;
    }
    
  std::ofstream fout(filename.toStdString().c_str());
  for(unsigned int i = 0; i < this->Sources.size(); ++i)
    {
    fout << "f " << this->Sources[0] << this->Sources[1] << std::endl; // 'f' stands for 'foreground'
    }
    
  for(unsigned int i = 0; i < this->Sinks.size(); ++i)
    {
    fout << "b " << this->Sinks[0] << this->Sinks[1] << std::endl; // 'b' stands for 'background'
    }
 
  fout.close();
}

void LidarSegmentationWidget::on_action_Selections_LoadForegroundSelectionsFromImage_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this,
     "Open Image", ".", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  typedef  itk::ImageFileReader< MaskImageType  > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.toStdString());
  reader->Update();

  std::vector<itk::Index<2> > pixels = Helpers::GetNonZeroPixels(reader->GetOutput());

  this->Sources = pixels;

  UpdateSelections();
}

void LidarSegmentationWidget::on_action_Selections_LoadBackgroundSelectionsFromImage_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this,
     "Open Image", ".", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  typedef  itk::ImageFileReader< MaskImageType  > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.toStdString());
  reader->Update();

  std::vector<itk::Index<2> > pixels = Helpers::GetNonZeroPixels(reader->GetOutput());

  this->Sinks = pixels;

  UpdateSelections();
}

void LidarSegmentationWidget::on_btnGenerateNeighborSinks_clicked()
{
  GenerateNeighborSinks();
}

void LidarSegmentationWidget::on_btnErodeSources_clicked()
{
  MaskImageType::Pointer sourcesImage = MaskImageType::New();
  sourcesImage->SetRegions(this->ImageRegion);
  Helpers::IndicesToBinaryImage(this->Sources, sourcesImage);

  typedef itk::BinaryBallStructuringElement<MaskImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElementBig;
  structuringElementBig.SetRadius(3);
  structuringElementBig.CreateStructuringElement();

  typedef itk::BinaryErodeImageFilter <MaskImageType, MaskImageType, StructuringElementType> BinaryErodeImageFilterType;

  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  erodeFilter->SetInput(sourcesImage);
  erodeFilter->SetKernel(structuringElementBig);
  erodeFilter->Update();

  //this->Sources.clear();

  this->Sources = Helpers::GetNonZeroPixels(erodeFilter->GetOutput());

  UpdateSelections();
}

void LidarSegmentationWidget::GenerateNeighborSinks()
{
  MaskImageType::Pointer sourcesImage = MaskImageType::New();
  sourcesImage->SetRegions(this->ImageRegion);
  Helpers::IndicesToBinaryImage(this->Sources, sourcesImage);
  Helpers::WriteImage(sourcesImage.GetPointer(), "sourcesImage.png");
  
  // Dilate the mask
  std::cout << "Dilating mask..." << std::endl;
  typedef itk::BinaryBallStructuringElement<MaskImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(1);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType> BinaryDilateImageFilterType;

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

  Helpers::WriteImage(xorFilter->GetOutput(), "boundaryOfSegmentation.png");
  
  // Iterate over the border pixels. If the closest pixel in the original segmentation has
  // a depth greater than a threshold, mark it as a new sink. Else, do not.
  std::cout << "Determining which boundary pixels should be declared background..." << std::endl;
  //std::cout << "There should be " << Helpers::CountNonZeroPixels(xorFilter->GetOutput())
  //          << " considered." << std::endl;
  typedef std::vector<itk::Index<2> > VectorOfPixelsType;
  VectorOfPixelsType newSinks;

  typedef itk::VectorIndexSelectionCastImageFilter<ImageType, FloatScalarImageType> IndexSelectionType;
  IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(3);
  indexSelectionFilter->SetInput(this->Image);
  indexSelectionFilter->Update();
  
  FloatScalarImageType::Pointer depthImage = indexSelectionFilter->GetOutput();

  //float sameObjectThreshold = 0.1f;

  
  VectorOfPixelsType consideredPixels;

  itk::ImageRegionIterator<MaskImageType> imageIterator(xorFilter->GetOutput(),
                                                        xorFilter->GetOutput()->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // If the current pixel is in question
      {
      consideredPixels.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }

  std::cout << "There are " << consideredPixels.size() << " potential new sink pixels." << std::endl;
  
  for(VectorOfPixelsType::const_iterator iter = consideredPixels.begin(); iter != consideredPixels.end(); ++iter)
    {
    //std::cout << "Considering pixel " << consideredCounter << " (index "
    //          << imageIterator.GetIndex() << ")" << std::endl;
    ImageType::PixelType currentPixel = this->Image->GetPixel(*iter);

    unsigned int radius = this->txtBackgroundCheckRadius->text().toUInt();
    ImageType::RegionType desiredRegion = Helpers::GetRegionInRadiusAroundPixel(*iter, radius);
    //std::cout << "desiredRegion: " << desiredRegion << std::endl;
    
    itk::ImageRegionIterator<MaskImageType> sourcesImageIterator(sourcesImage, desiredRegion);
    
    std::vector<float> nonForegroundDepths;
    std::vector<float> foregroundDepths;
    while(!sourcesImageIterator.IsAtEnd())
      {
      if(sourcesImageIterator.Get())
        {
        foregroundDepths.push_back(depthImage->GetPixel(sourcesImageIterator.GetIndex()));
        }
      else
        {
        nonForegroundDepths.push_back(depthImage->GetPixel(sourcesImageIterator.GetIndex()));
        }
      ++sourcesImageIterator;
      }

    float nonForegroundMedian = Helpers::VectorMedian(nonForegroundDepths);
    float foregroundMedian = Helpers::VectorMedian(foregroundDepths);

    float difference = fabs(foregroundMedian - nonForegroundMedian);

    //if(difference > this->txtBackgroundThreshold->text().toFloat())
    if(difference > this->txtBackgroundThreshold->text().toFloat())
      {
      std::cout << "Difference was " << difference << " so this is a sink pixel." << std::endl;
      newSinks.push_back(*iter);
      }
    else
      {
      std::cout << "Difference was " << difference << " so this is NOT a sink pixel." << std::endl;
      }
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
  Helpers::WriteImage(newSinksImage.GetPointer(), "newSinks.png");

  //std::cout << "Out of " << consideredCounter << " pixels considered, "
  //          << backgroundCounter << " were declared background." << std::endl;
  // Set the new sinks
  std::cout << "Setting " << newSinks.size() << " new sinks." << std::endl;

  // Modify the list of sinks so it can be retrieved by the MainWindow after the segmentation is finished
  this->Sinks.insert(this->Sinks.end(), newSinks.begin(), newSinks.end());

  UpdateSelections();
}

void LidarSegmentationWidget::on_action_Selections_LoadSelectionsFromImage_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this,
     "Open Image", ".", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  typedef  itk::ImageFileReader< RGBImageType  > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.toStdString());
  reader->Update();

  RGBPixelType greenPixel;
  greenPixel.SetRed(0);
  greenPixel.SetGreen(255);
  greenPixel.SetBlue(0);
  
  RGBPixelType redPixel;
  redPixel.SetRed(255);
  redPixel.SetGreen(0);
  redPixel.SetBlue(0);
  
  itk::ImageRegionConstIterator<RGBImageType> imageIterator(reader->GetOutput(),
                                                            reader->GetOutput()->GetLargestPossibleRegion());
 
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == greenPixel)
      {
      this->Sources.push_back(imageIterator.GetIndex());
      }
    else if(imageIterator.Get() == redPixel)
      {
      this->Sinks.push_back(imageIterator.GetIndex());
      }
 
    ++imageIterator;
    }

  UpdateSelections();
}


void LidarSegmentationWidget::on_action_Selections_LoadSelectionsFromText_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this,
     "Open Image", ".", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }
    
  std::ifstream fin(filename.toStdString().c_str());
 
  if(fin == NULL)
    {
    std::cout << "Cannot open file." << std::endl;
    }
 
  std::string line;
 
  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    char selectionType;
    itk::Index<2> pixel;
    ss >> selectionType >> pixel[0] >> pixel[1];
    if(selectionType == 'f')
      {
      this->Sources.push_back(pixel);
      }
    else if(selectionType == 'b')
      {
      this->Sinks.push_back(pixel);
      }
    else
      {
      std::stringstream ss;
      ss << "Selectiontype is " << selectionType << " - should be 'f' or 'b' (foreground or background)";
      throw std::runtime_error(ss.str());
      }
    }
}


void LidarSegmentationWidget::on_btnCut_clicked()
{
  // Normalize the image
  ImageType::Pointer normalizedImage = ImageType::New();
  normalizedImage->SetNumberOfComponentsPerPixel(this->Image->GetNumberOfComponentsPerPixel());
  normalizedImage->SetRegions(this->Image->GetLargestPossibleRegion());
  normalizedImage->Allocate();

  std::cout << "Normalizing image..." << std::endl;
  Helpers::NormalizeImage(this->Image.GetPointer(), normalizedImage.GetPointer());

  std::cout << "Normalized image has " << normalizedImage->GetNumberOfComponentsPerPixel() << " channels." << std::endl;
  
  //this->GraphCut.SetImage(this->Image);
  this->GraphCut.SetImage(normalizedImage.GetPointer());

  this->GraphCut.Debug = this->chkDebug->isChecked();
  //this->GraphCut.SecondStep = this->chkSecondStep->isChecked();
  
  this->GraphCut.IncludeDepthInHistogram = this->chkDepthHistogram->isChecked();
  this->GraphCut.IncludeColorInHistogram = this->chkColorHistogram->isChecked();

  this->GraphCut.BackgroundThreshold = this->txtBackgroundThreshold->text().toDouble();
  
  if(this->GraphCut.DifferenceFunction)
    {
    delete this->GraphCut.DifferenceFunction;
    }

  // Setup the Difference object
  if(this->chkDepthDifference->isChecked() && !this->chkColorDifference->isChecked())
    {
    std::cout << "Using depth-only N-weights." << std::endl;
    this->GraphCut.DifferenceFunction = new DepthDifference;
    }
  else if(!this->chkDepthDifference->isChecked() && this->chkColorDifference->isChecked())
    {
    std::cout << "Using color-only N-weights." << std::endl;
    this->GraphCut.DifferenceFunction = new ColorDifference;
    }
  else if(this->chkDepthDifference->isChecked() && this->chkColorDifference->isChecked())
    {
    std::cout << "Using depth+color N-weights." << std::endl;
    std::vector<float> weights(4,1.0f);
    weights[0] = spinRWeight->value();
    weights[1] = spinGWeight->value();
    weights[2] = spinBWeight->value();
    weights[3] = spinDWeight->value();

    std::cout << "Weights: " << weights[0] << " " << weights[1] << " " << weights[2] << " " << weights[3] << std::endl;

    this->GraphCut.DifferenceFunction = new WeightedDifference(weights);
    }
  else
    {
    std::stringstream ss;
    ss << "Something is wrong - you must select depth, color, or both." << std::endl;
    throw std::runtime_error(ss.str());
    }

  // Get the number of bins from the slider
  this->GraphCut.SetNumberOfHistogramBins(this->sldHistogramBins->value());

  if(this->sldLambda->value() == 0)
    {
    QMessageBox msgBox;
    msgBox.setText("You must select lambda > 0!");
    msgBox.exec();
    return;
    }

  // Setup the graph cut from the GUI and the scribble selection
  this->GraphCut.SetLambda(ComputeLambda());
  this->GraphCut.SetSources(this->Sources);
  this->GraphCut.SetSinks(this->Sinks);
  //this->GraphCut.SetSources(this->LeftInteractorStyle->GetForegroundSelection());
  //this->GraphCut.SetSinks(this->LeftInteractorStyle->GetBackgroundSelection());

  // Setup and start the actual cut computation in a different thread
  QFuture<void> future = QtConcurrent::run(this->GraphCut, &ImageGraphCut::PerformSegmentation);
  this->FutureWatcher.setFuture(future);

  this->ProgressDialog->exec();

}

#if 0
void InnerWidget::actionFlip_Image_triggered()
{
  this->CameraUp[1] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->Refresh();
}
#endif

#if 0
void InnerWidget::actionSave_Segmentation_triggered()
{
  // Ask the user for a filename to save the segment mask image to

  QString fileName = QFileDialog::getSaveFileName(this,
    tr("Save Segment Mask Image"), "/home/doriad", tr("Image Files (*.png *.bmp)"));
/*
  // Convert the image from a 1D vector image to an unsigned char image
  typedef itk::CastImageFilter< GrayscaleImageType,
                                itk::Image<itk::CovariantVector<unsigned char, 1>, 2 > > CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(this->GraphCut->GetSegmentMask());

  typedef itk::NthElementImageAdaptor< itk::Image<itk:: CovariantVector<unsigned char, 1>, 2 >,
    unsigned char> ImageAdaptorType;

  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SelectNthElement(0);
  adaptor->SetImage(castFilter->GetOutput());
*/
  // Write the file
  //typedef  itk::ImageFileWriter< ImageAdaptorType > WriterType;
  typedef  itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  //writer->SetInput(adaptor);
  writer->SetInput(this->GraphCut->GetSegmentMask());
  writer->Update();

}
#endif

void LidarSegmentationWidget::OpenFile(const std::string& fileName)
{
  // Store the path so that temp files can be written to the same path as the input image
  QFileInfo fileInfo(fileName.c_str());
  std::string workingDirectory = fileInfo.absoluteDir().absolutePath().toStdString() + "/";
  std::cout << "Working directory set to: " << workingDirectory << std::endl;
  QDir::setCurrent(QString(workingDirectory.c_str()));
  
  // Read file
  itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();

  reader->SetFileName(fileName);
  reader->Update();

  if(reader->GetOutput()->GetNumberOfComponentsPerPixel() < 4)
    {
    std::cerr << "The input image has " << reader->GetOutput()->GetNumberOfComponentsPerPixel()
              << " components, but (at least) 4 are required." << std::endl;
    return;
    }

  Helpers::DeepCopy(reader->GetOutput(), this->Image.GetPointer());

  // Store the region so we can access it without needing to care which image it comes from
  this->ImageRegion = reader->GetOutput()->GetLargestPossibleRegion();

  // Clear everything
  //this->LeftRenderer->RemoveAllViewProps();
  //this->RightRenderer->RemoveAllViewProps();
  this->Sources.clear();
  this->Sinks.clear();

  UpdateSelections();

  // Convert the ITK image to a VTK image and display it
  Helpers::ITKImageToVTKImage(reader->GetOutput(), this->OriginalImageData);

  this->OriginalImageSliceMapper->SetInputData(this->OriginalImageData);
  this->OriginalImageSlice->SetMapper(this->OriginalImageSliceMapper);

  //this->LeftRenderer->AddViewProp(this->OriginalImageSlice);
  this->LeftRenderer->ResetCamera();

  // Setup the scribble canvas
  Helpers::SetImageSize(this->OriginalImageData, this->SourceSinkImageData);
  Helpers::CreateTransparentImage(this->SourceSinkImageData);

  this->LeftSourceSinkImageSliceMapper->SetInputData(this->SourceSinkImageData);
  this->LeftSourceSinkImageSlice->SetMapper(this->LeftSourceSinkImageSliceMapper);

  this->RightSourceSinkImageSliceMapper->SetInputData(this->SourceSinkImageData);
  this->RightSourceSinkImageSlice->SetMapper(this->RightSourceSinkImageSliceMapper);

  // If this is called, the image disappears in the *left* renderer???
  this->RightRenderer->AddViewProp(this->RightSourceSinkImageSlice); 
  this->RightRenderer->ResetCamera();

  // This also adds the ImageSlice to the renderers
  this->LeftInteractorStyle->InitializeTracer(this->LeftSourceSinkImageSlice);

  // This is done inside the InitializeTracer call
  //this->LeftRenderer->AddViewProp(this->LeftSourceSinkImageSlice); 

  this->Refresh();
}

void LidarSegmentationWidget::Refresh()
{
  std::cout << "Refresh()" << std::endl;
  //this->SourceSinkImageSliceMapper->Render();
  //this->SourceSinkImageSliceMapper->Modified();
  
  this->qvtkWidgetRight->GetRenderWindow()->Render();
  this->qvtkWidgetLeft->GetRenderWindow()->Render();
  this->qvtkWidgetRight->GetInteractor()->Render();
  this->qvtkWidgetLeft->GetInteractor()->Render();
}

void LidarSegmentationWidget::on_action_View_DepthImage_triggered()
{
  //Helpers::ITKImageChannelToVTKImage(thresholdFilter->GetOutput(), 3, this->OriginalImageData);
  
  typedef itk::Image<ImageType::InternalPixelType, 2> ScalarImageType;
  ScalarImageType::Pointer depthImage = ScalarImageType::New();
  
  Helpers::ExtractChannel(this->Image.GetPointer(), 3, depthImage.GetPointer());
  
  typedef itk::ThresholdImageFilter<ScalarImageType> ThresholdImageFilterType;

  float maxDepth = 20.0f;
  ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
  thresholdFilter->SetInput(depthImage);
  thresholdFilter->ThresholdAbove(maxDepth);
  thresholdFilter->SetOutsideValue(maxDepth);
  thresholdFilter->Update();
  
  Helpers::ITKScalarImageToVTKImage(thresholdFilter->GetOutput(), this->OriginalImageData);
  
  Refresh();
}

void LidarSegmentationWidget::on_action_View_ColorImage_triggered()
{
  Helpers::ITKImageToVTKImage(this->Image.GetPointer(), this->OriginalImageData);
  Refresh();
}

void LidarSegmentationWidget::on_action_Export_InputScreenshot_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this, "Save PNG Screenshot", "", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(this->qvtkWidgetLeft->GetRenderWindow());
  //windowToImageFilter->SetMagnification(3);
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(filename.toStdString().c_str());
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
}

void LidarSegmentationWidget::on_action_Export_ResultScreenshot_triggered()
{
  QString filename = QFileDialog::getSaveFileName(this, "Save PNG Screenshot", "", "PNG Files (*.png)");

  if(filename.isEmpty())
    {
    return;
    }

  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(this->qvtkWidgetRight->GetRenderWindow());
  //windowToImageFilter->SetMagnification(3);
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(filename.toStdString().c_str());
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
}

void LidarSegmentationWidget::on_btnReseedForeground_clicked()
{
  std::vector<itk::Index<2> > pixels = Helpers::GetNonZeroPixels(this->GraphCut.GetSegmentMask());

  this->Sources = pixels;

  UpdateSelections();
}

void LidarSegmentationWidget::on_btnShowStrokesRight_clicked()
{
  RightSourceSinkImageSlice->VisibilityOn();
  Refresh();
}

void LidarSegmentationWidget::on_btnHideStrokesRight_clicked()
{
  RightSourceSinkImageSlice->VisibilityOff();
  Refresh();
}


void LidarSegmentationWidget::on_btnShowStrokesLeft_clicked()
{
  LeftSourceSinkImageSlice->VisibilityOn();
  Refresh();
}

void LidarSegmentationWidget::on_btnHideStrokesLeft_clicked()
{
  LeftSourceSinkImageSlice->VisibilityOff();
  Refresh();
}
