#ifndef VOITOMASK_H
#define VOITOMASK_H

#include "itkImage.h"
#include "VoiFileIO.h"

namespace nih {

template<typename PixelType>
bool ConvertVoiToMask2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, const VoiFile::Slice &stSlice, bool bErase = true);

template<typename PixelType>
bool ConvertVoiToMask3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, const VoiFile &clVoi, bool bErase = true);

} // end namespace nih

#include "VoiToMask.hpp"

#endif // !VOITOMASK_H
