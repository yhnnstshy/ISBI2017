#SET(CMAKE_PREFIX_PATH "/home/tsehayyk/Work/ProstateCAD/Source/Dependencies/Linux_GCC4.9/Release/InsightToolkit-4.7.1/lib/cmake/ITK-4.7")

SET(CMAKE_CXX_FLAGS "-std=c++11")

FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKGDCM ITKIOGDCM ITKIOMeta ITKIONIFTI ITKIONRRD ITKIOJPEG ITKIOPNG ITKMesh ITKQuadEdgeMesh)

IF (NOT ITK_FOUND)
  MESSAGE(FATAL "ITK is required.")
ENDIF()

INCLUDE(${ITK_USE_FILE})

ADD_LIBRARY(CADLib
  CAD/Common.h CAD/Common.cpp
  CAD/CodeBookTree.h CAD/CodeBookTree.hpp
  CAD/VoiFileIO.h CAD/VoiFileIO.cpp
  CAD/VoiToMask.h CAD/VoiToMask.hpp
  CAD/MedicalImages2HED.h CAD/MedicalImages2HED.cpp
)


add_executable(runhed CAD/RunHED.cpp CAD/bsdgetopt.h CAD/bsdgetopt.c)
target_link_libraries(runhed CADLib ${Caffe_LINK} ${ITK_LIBRARIES} "-lhdf5_cpp -lhdf5")
  

caffe_default_properties(runhed)

# set back RUNTIME_OUTPUT_DIRECTORY
set_target_properties(runhed PROPERTIES
RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/examples/CAD")
caffe_set_solution_folder(runhed examples)


# install
install(TARGETS runhed DESTINATION bin)


if(UNIX OR APPLE)
  # Funny command to make tutorials work
  # TODO: remove in future as soon as naming is standartaized everywhere
  set(__outname ${PROJECT_BINARY_DIR}/examples/CAD/runhed${Caffe_POSTFIX})
  add_custom_command(TARGET runhed POST_BUILD
                     COMMAND ln -sf "${__outname}" "${__outname}.bin")
endif()

