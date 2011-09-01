/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

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

/* This is the main GUI class of this project. It is a QMainWindow
 * so that we can use a File menu. It contains an instance of our main functional
 * class ImageGraphCutBase and our custom scribble interactor style vtkGraphCutInteractorStyle.
 * It also contains a CProgressThread so that we can display a progress bar in marquee
 * mode during long computations.
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// Qt
#include "ui_MainWindow.h"

// Custom
#include "ProgressThread.h"
class InteractorStyleScribble;
class InteractorStyleImageNoLevel;

// Forward declarations
class vtkImageSlice;
class vtkImageSliceMapper;

class MainWindow : public QMainWindow, private Ui::MainWindow
{
Q_OBJECT
public:
  MainWindow(QWidget *parent = 0);

public slots:
  // Menu items
  void on_actionOpenImage_triggered();
  void on_actionSaveSegmentation_triggered();
  void on_actionFlipImage_triggered();
  void on_actionExit_triggered();
  void on_actionSaveSelectionsAsImage_triggered();
  void on_actionSaveSelectionsAsText_triggered();
  void on_actionLoadSelectionsFromImage_triggered();
  void on_actionLoadSelectionsFromText_triggered();

  // Buttons, radio buttons, and sliders
  void on_btnClearSelections_clicked();
  void on_btnClearBackground_clicked();
  void on_btnClearForeground_clicked();
  void on_btnCut_clicked();
  void on_radForeground_clicked();
  void on_radBackground_clicked();
  void sldHistogramBins_valueChanged();
  
  // Setting lambda must be handled specially because we need to multiply the percentage set by the slider by the MaxLambda set in the text box
  void UpdateLambda();

  // These slots handle running the progress bar while the computations are done in a separate thread.
  void StartProgressSlot();
  void StopProgressSlot();

  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.
  void OpenFile();
  
  void UpdateSelections();
  
protected:

  void DisplaySegmentationResult();
  
  // A class to do the main computations in a separate thread so we can display a marquee progress bar.
  ProgressThread<ImageType> SegmentationThread;

  // Compute lambda by multiplying the percentage set by the slider by the MaxLambda set in the text box.
  float ComputeLambda();

  // Left pane
  vtkSmartPointer<InteractorStyleScribble> LeftInteractorStyle;
  vtkSmartPointer<vtkImageSliceMapper> OriginalImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> OriginalImageSlice;
  vtkSmartPointer<vtkRenderer> LeftRenderer;
  vtkSmartPointer<vtkImageData> OriginalImageData;
  
  // Right pane
  vtkSmartPointer<vtkImageData> ResultImageData;
  vtkSmartPointer<vtkImageSliceMapper> ResultImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> ResultImageSlice;
  vtkSmartPointer<vtkRenderer> RightRenderer;
  vtkSmartPointer<InteractorStyleImageNoLevel> RightInteractorStyle;
  
  // Both panes
  vtkSmartPointer<vtkImageData> SourceSinkImageData; // This data can be used by both the Left and Right SourceSinkImageSliceMapper
  
  vtkSmartPointer<vtkImageSliceMapper> LeftSourceSinkImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> LeftSourceSinkImageSlice;
  
  vtkSmartPointer<vtkImageSliceMapper> RightSourceSinkImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> RightSourceSinkImageSlice;
  
  // Refresh both renderers and render windows
  void Refresh();

  // The main segmentation class. This will be instantiated as a ImageGraphCut after the user selects whether to open a color or grayscale image.
  ImageGraphCut GraphCut;

  // Allows the background color to be changed
  double BackgroundColor[3];

  // We set this when the image is opeend. We sometimes need to know how big the image is.
  itk::ImageRegion<2> ImageRegion;
  
  void ScribbleEventHandler(vtkObject* caller, long unsigned int eventId, void* callData);
  
  std::vector<itk::Index<2> > Sources;
  std::vector<itk::Index<2> > Sinks;
  
  bool Flipped;
  void SetCameraPosition1();
  void SetCameraPosition2();
};

#endif
