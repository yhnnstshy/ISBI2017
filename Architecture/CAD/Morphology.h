#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include "itkImage.h"

namespace nih {

template<typename PixelType, typename StructuringElementType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer Dilate(typename itk::Image<PixelType, Dimension>::Pointer p_clImage, const StructuringElementType &clKernel);

template<typename PixelType, typename StructuringElementType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer Erode(typename itk::Image<PixelType, Dimension>::Pointer p_clImage, const StructuringElementType &clKernel);

// Box-based dilation
template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer Dilate2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY);

// Cube-based dilation
template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer Dilate3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY, unsigned int uiRadiusZ);

// Box-based erode
template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer Erode2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY);

// Cube-based dilation
template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer Erode3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, unsigned int uiRadiusX, unsigned int uiRadiusY, unsigned int uiRadiusZ);

} // end namespace nih

#include "Morphology.hpp"

#endif // MORPHOLOGY_H
