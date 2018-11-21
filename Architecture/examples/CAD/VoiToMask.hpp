#include "VoiToMask.h"

namespace nih {

inline bool HorizontalRayIntersects(float x, float y, const VoiFile::PointType &a, const VoiFile::PointType &b) {
  const float fSmall = 1e-10f;
  const float fDiffY = a[1] - b[1];

  if (std::abs(fDiffY) < fSmall) {
    if (std::abs(a[1] - y) < fSmall && a[0] >= x)
      return true;

    // NOTE: Contour is supposedly close, b[1] will be a[1] in a different contour
    //if (std::abs(b[1] - y) < fSmall && b[0] >= x)
      //return true;
  }
  else {
    const float t = (y - b[1])/fDiffY;
    const float fXInt = a[0]*t + b[0]*(1.0f - t);

    if (t > 0.0f && t <= 1.0f && fXInt >= x)
      return true;
  }

  return false;
}

template<typename PixelType>
bool ConvertVoiToMask2D(typename itk::Image<PixelType, 2>::Pointer p_clImage, const VoiFile::Slice &stSlice, bool bErase) {
  typedef itk::Image<PixelType, 2> ImageType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  //typedef typename ImageType::PointType PointType;
  typedef typename ImageType::IndexType IndexType;

  if (!p_clImage)
    return false;

  const SizeType &clSize = p_clImage->GetBufferedRegion().GetSize();
  //const SpacingType &clSpacing = p_clImage->GetSpacing();
  //const PointType &clOrigin = p_clImage->GetOrigin();

  if (clSize[0] <= 0 || clSize[1] <= 0 || clSize[2] <= 0)
    return false;

  typedef VoiFile::PointType VoiPointType;

  const std::vector<std::vector<VoiPointType> > &vContours = stSlice.vContours;

  for (size_t i = 0; i < vContours.size(); ++i) {
    const std::vector<VoiPointType> &vContour = vContours[i];

    for (int y = 0; y < clSize[1]; ++y) {
      //const float fY = (float)(clSpacing[1]*y);
      const float fY = (float)y;

      for (int x = 0; x < clSize[0]; ++x) {
        //const float fX = (float)(clSpacing[0]*x);
        const float fX = (float)x;

        const IndexType clIndex = { x, y };

        unsigned int uiIntersectCount = 0;

        // Count the number of intersections
        for (size_t j = 1; j < vContour.size(); ++j) {
          if (HorizontalRayIntersects(fX, fY, vContour[j-1], vContour[j]))
            ++uiIntersectCount;
        }

        if (HorizontalRayIntersects(fX, fY, vContour.back(), vContour.front()))
          ++uiIntersectCount;

        if (uiIntersectCount & 1)
          p_clImage->SetPixel(clIndex, (PixelType)1);
        else if (bErase)
          p_clImage->SetPixel(clIndex, (PixelType)0);
      }
    }


  }

  return true;
}

template<typename PixelType>
bool ConvertVoiToMask3D(typename itk::Image<PixelType, 3>::Pointer p_clImage, const VoiFile &clVoi, bool bErase) {
  typedef itk::Image<PixelType, 3> ImageType3D;
  typedef itk::Image<PixelType, 2> ImageType2D;

  typedef typename ImageType3D::SizeType SizeType3D;
  typedef typename ImageType2D::SizeType SizeType2D;

  typedef typename ImageType3D::SpacingType SpacingType3D;
  typedef typename ImageType2D::SpacingType SpacingType2D;

  if (!p_clImage)
    return false;

  const SizeType3D &clSize = p_clImage->GetBufferedRegion().GetSize();
  const SpacingType3D &clSpacing = p_clImage->GetSpacing();

  typedef VoiFile::PointType VoiPointType;

  const std::vector<VoiFile::Slice> &vSlices = clVoi.GetSlices();

  PixelType * const p_buffer = p_clImage->GetPixelContainer()->GetImportPointer();

  const typename SizeType2D::SizeValueType sliceSize = clSize[0]*clSize[1];

  int iSucceed = 1;

#pragma omp parallel
  {
    typename ImageType2D::Pointer p_clSlice = ImageType2D::New();
   
    {
      const SizeType2D clSize2D = { clSize[0], clSize[1] };
      SpacingType2D clSpacing2D;
   
      clSpacing2D[0] = clSpacing[0];
      clSpacing2D[1] = clSpacing[1];
   
      p_clSlice->SetRegions(clSize2D);
      p_clSlice->SetSpacing(clSpacing2D);
    }

    bool bThreadSucceed = true;
    
#pragma omp for schedule(dynamic) reduction(&:iSucceed)
    for (int i = 0; i < (int)vSlices.size(); ++i) {
      // XXX: Uhh, not sure if this starts at 0 or 1
      const VoiFile::Slice &stSlice = vSlices[i];
      const unsigned int uiSliceNumber = stSlice.uiSliceNumber;
  
      if (uiSliceNumber >= clSize[2]) {
        iSucceed = 0;
        continue;
      }
  
      p_clSlice->GetPixelContainer()->SetImportPointer(p_buffer + uiSliceNumber * sliceSize, sliceSize, false);
  
      ConvertVoiToMask2D<PixelType>(p_clSlice, stSlice, bErase);
    }

  } // end parallel section

  return iSucceed;
}

} // end namespace nih
