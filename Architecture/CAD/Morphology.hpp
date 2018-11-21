#include "Morphology.h"
#include "itkFlatStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"

namespace nih {

template<typename PixelType, typename StructuringElementType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer Dilate(typename itk::Image<PixelType, Dimension>::Pointer p_clImage, const StructuringElementType &clKernel) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType, StructuringElementType> DilateFilterType;

  typename DilateFilterType::Pointer p_clFilter = DilateFilterType::New();

  if (!p_clFilter)
    return typename ImageType::Pointer();

  p_clFilter->SetInput(p_clImage);
  p_clFilter->SetKernel(clKernel);
  
  try {
    p_clFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  typename ImageType::Pointer p_clOutputImage = p_clFilter->GetOutput();

  return p_clOutputImage;
}

template<typename PixelType, typename StructuringElementType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer Erode(typename itk::Image<PixelType, Dimension>::Pointer p_clImage, const StructuringElementType &clKernel) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, StructuringElementType> DilateFilterType;

  typename DilateFilterType::Pointer p_clFilter = DilateFilterType::New();

  if (!p_clFilter)
    return typename ImageType::Pointer();

  p_clFilter->SetInput(p_clImage);
  p_clFilter->SetKernel(clKernel);
  
  try {
    p_clFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  typename ImageType::Pointer p_clOutputImage = p_clFilter->GetOutput();

  return p_clOutputImage;
}

template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer Dilate2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY) {
  enum { Dimension = 2 };
  typedef itk::FlatStructuringElement<Dimension> StructuringElementType;
  typedef typename StructuringElementType::RadiusType RadiusType;

  RadiusType clRadius;

  clRadius[0] = uiRadiusX;
  clRadius[1] = uiRadiusY;
  
  StructuringElementType clBox = StructuringElementType::Box(clRadius);

  return Dilate<PixelType, StructuringElementType, Dimension>(p_clImage, clBox);
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer Dilate3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY, unsigned int uiRadiusZ) {
  enum { Dimension = 3 };
  typedef itk::FlatStructuringElement<Dimension> StructuringElementType;
  typedef typename StructuringElementType::RadiusType RadiusType;

  RadiusType clRadius;

  clRadius[0] = uiRadiusX;
  clRadius[1] = uiRadiusY;
  clRadius[2] = uiRadiusZ;
  
  StructuringElementType clBox = StructuringElementType::Box(clRadius);

  return Dilate<PixelType, StructuringElementType, Dimension>(p_clImage, clBox);
}

template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer Erode2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY) {
  enum { Dimension = 2 };
  typedef itk::FlatStructuringElement<Dimension> StructuringElementType;
  typedef typename StructuringElementType::RadiusType RadiusType;

  RadiusType clRadius;

  clRadius[0] = uiRadiusX;
  clRadius[1] = uiRadiusY;
  
  StructuringElementType clBox = StructuringElementType::Box(clRadius);

  return Erode<PixelType, StructuringElementType, Dimension>(p_clImage, clBox);
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer Erode3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY, unsigned int uiRadiusZ) {
  enum { Dimension = 3 };
  typedef itk::FlatStructuringElement<Dimension> StructuringElementType;
  typedef typename StructuringElementType::RadiusType RadiusType;

  RadiusType clRadius;

  clRadius[0] = uiRadiusX;
  clRadius[1] = uiRadiusY;
  clRadius[2] = uiRadiusZ;
  
  StructuringElementType clBox = StructuringElementType::Box(clRadius);

  return Erode<PixelType, StructuringElementType, Dimension>(p_clImage, clBox);
}

} // end namespace nih
