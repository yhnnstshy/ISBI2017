SET(CMAKE_CXX_FLAGS "-std=c++11")

FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKGDCM ITKIOGDCM ITKIOMeta ITKIONIFTI ITKIONRRD ITKIOJPEG ITKIOPNG ITKMesh ITKQuadEdgeMesh)

IF (NOT ITK_FOUND)
  MESSAGE(FATAL "ITK is required.")
ENDIF()

INCLUDE(${ITK_USE_FILE})

ADD_LIBRARY(CADLib
  Common.h Common.cpp
  ProstateCAD.h ProstateCAD.cpp
  ROCCurve.h    ROCCurve.cpp
  CodeBookTree.h CodeBookTree.hpp
  VoiFileIO.h VoiFileIO.cpp
  VoiToMask.h VoiToMask.hpp
  MedicalImages2HED.h MedicalImages2HED.cpp
  PiecewiseLinearInterpolation.h
  Morphology.h Morphology.hpp
  MeanShift.h MeanShift.hpp
)


add_executable(Detect Detect.cpp bsdgetopt.h bsdgetopt.c)
target_link_libraries(Detect CADLib ${Caffe_LINK} ${ITK_LIBRARIES} "-lhdf5_cpp -lhdf5")
  
ADD_EXECUTABLE(EvalFROC EvalGleasonFROC.cpp bsdgetopt.h bsdgetopt.c)
TARGET_LINK_LIBRARIES(EvalFROC CADLib ${ITK_LIBRARIES})

ADD_EXECUTABLE(EvalROC EvalProstateCAD.cpp bsdgetopt.h bsdgetopt.c)
TARGET_LINK_LIBRARIES(EvalROC CADLib ${ITK_LIBRARIES})

ADD_EXECUTABLE(AverageFROCCurves AverageFROCCurves.cpp bsdgetopt.h bsdgetopt.c)
TARGET_LINK_LIBRARIES(AverageFROCCurves CADLib ${ITK_LIBRARIES})

ADD_EXECUTABLE(AverageROCCurves AverageROCCurves.cpp bsdgetopt.h bsdgetopt.c)
TARGET_LINK_LIBRARIES(AverageROCCurves CADLib ${ITK_LIBRARIES})

caffe_default_properties(Detect)

# set back RUNTIME_OUTPUT_DIRECTORY
set_target_properties(Detect PROPERTIES
RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/CAD")
caffe_set_solution_folder(Detect CAD)
