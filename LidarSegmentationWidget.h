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

#ifndef LidarSegmentationWidget_H
#define LidarSegmentationWidget_H

// GUI
#include "ui_LidarSegmentationWidget.h"

// Qt
#include <QFutureWatcher>
#include <QProgressDialog>

// VTK
#include <vtkSmartPointer.h>

// Custom
class InteractorStyleScribble;
class InteractorStyleImageNoLevel;
#include "ImageGraphCut.h"

// Forward declarations
class vtkImageSlice;
class vtkImageSliceMapper;

class LidarSegmentationWidget : public QMainWindow, private Ui::LidarSegmentationWidget
{
Q_OBJECT
public:
  LidarSegmentationWidget(QWidget *parent = 0);
  LidarSegmentationWidget(const std::string& fileName);
  
public slots:
  // Menu items
  void on_actionOpenImage_triggered();
  void on_actionSaveSegmentation_triggered();
  void on_actionFlipImage_triggered();
  void on_actionExit_triggered();

  // Selections menu
  void on_action_Selections_SaveSelectionsAsImage_triggered();
  void on_action_Selections_SaveForegroundSelectionsAsImage_triggered();
  void on_action_Selections_SaveBackgroundSelectionsAsImage_triggered();
  void on_action_Selections_SaveSelectionsAsText_triggered();
  void on_action_Selections_LoadSelectionsFromImage_triggered();
  void on_action_Selections_LoadSelectionsFromText_triggered();
  void on_action_Selections_LoadForegroundSelectionsFromImage_triggered();
  void on_action_Selections_LoadBackgroundSelectionsFromImage_triggered();

  // View menu
  void on_action_View_DepthImage_triggered();
  void on_action_View_ColorImage_triggered();

  // Export menu
  void on_action_Export_InputScreenshot_triggered();
  void on_action_Export_ResultScreenshot_triggered();
  
  // Buttons, radio buttons, and sliders
  void on_btnGenerateNeighborSinks_clicked();
  void on_btnErodeSources_clicked();
  
  void on_btnClearSelections_clicked();
  void on_btnClearBackground_clicked();
  void on_btnClearForeground_clicked();
  void on_btnCut_clicked();
  void on_radForeground_clicked();
  void on_radBackground_clicked();
  void on_sldHistogramBins_valueChanged();
  
  /** Setting lambda must be handled specially because we need to multiply the percentage
   *  set by the slider by the MaxLambda set in the text box */
  void UpdateLambda();

  /** Perform an action when the segmentation has finished. */
  void slot_SegmentationComplete();

  /** Open the specified file as a greyscale or color image, depending on which type the user
   * has specified through the file menu.
   */
  void OpenFile(const std::string& fileName);
  
  void UpdateSelections();
  
private:
  /** A constructor that can be used by all other constructors. */
  void SharedConstructor();

  void GenerateNeighborSinks();
  
  void DisplaySegmentationResult();

  /** Compute lambda by multiplying the percentage set by the slider by the MaxLambda set in the text box. */
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
  
  /** Both panes - This data can be used by both the Left and Right SourceSinkImageSliceMapper */
  vtkSmartPointer<vtkImageData> SourceSinkImageData; 
  
  vtkSmartPointer<vtkImageSliceMapper> LeftSourceSinkImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> LeftSourceSinkImageSlice;
  
  vtkSmartPointer<vtkImageSliceMapper> RightSourceSinkImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> RightSourceSinkImageSlice;
  
  /** Refresh both renderers and render windows */
  void Refresh();

  /** The main segmentation class. This will be instantiated as a ImageGraphCut
   *  after the user selects an image. */
  ImageGraphCut GraphCut;

  /** Allows the background color to be changed */
  double BackgroundColor[3];

  /** We set this when the image is opeend. We sometimes need to know how big the image is. */
  itk::ImageRegion<2> ImageRegion;
  
  void ScribbleEventHandler(vtkObject* caller, long unsigned int eventId, void* callData);
  
  std::vector<itk::Index<2> > Sources;
  std::vector<itk::Index<2> > Sinks;
  
  bool Flipped;
  void SetCameraPosition1();
  void SetCameraPosition2();

  ImageType::Pointer Image;
  
  bool Debug;

  /** These members are to handle the progress bar during the segmentation.
   *  They have to be members so we can setup the connections once during the constructor.
   */
  QFutureWatcher<void> FutureWatcher;
  QProgressDialog* ProgressDialog;
};

#endif
