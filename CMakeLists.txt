cmake_minimum_required(VERSION 2.6)

Project(InteractiveLidarSegmentation)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Where to copy executables when 'make install' is run
SET( INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

QT4_WRAP_UI(LidarSegmentationUISrcs LidarSegmentationWidget.ui)
QT4_WRAP_CPP(LidarSegmentationMOCSrcs LidarSegmentationWidget.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${USE_VTK_FILE} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${USE_ITK_FILE} )

add_library(MaxFlow graph.cpp maxflow.cpp)

ADD_EXECUTABLE(InteractiveLidarSegmentation InteractiveLidarSegmentation.cpp
LidarSegmentationWidget.cpp
Helpers.cpp
Difference.cpp
InteractorStyleScribble.cxx
InteractorStyleImageNoLevel.cxx
ImageGraphCut.cxx
${LidarSegmentationMOCSrcs} ${LidarSegmentationUISrcs})
TARGET_LINK_LIBRARIES(InteractiveLidarSegmentation ${VTK_LIBRARIES}
QVTK ${QT_LIBRARIES} QtXml MaxFlow
${ITK_LIBRARIES}
)
INSTALL( TARGETS InteractiveLidarSegmentation RUNTIME DESTINATION ${INSTALL_DIR} )

# ADD_EXECUTABLE(NonInteractive NonInteractive.cpp
# graph.cpp maxflow.cpp
# Helpers.cpp
# )
# TARGET_LINK_LIBRARIES(NonInteractive vtkHybrid vtkWidgets
# ITKCommon ITKIO ITKStatistics
# )
