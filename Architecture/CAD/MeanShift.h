#ifndef MEANSHIFT_H
#define MEANSHIFT_H

#include <utility>
#include <vector>
#include "Common.h"

namespace nih {

template<typename PixelType, unsigned int Dimension>
class MeanShift {
public:
  enum { DIMENSION = Dimension };

  typedef itk::Image<PixelType, DIMENSION> ImageType;
  typedef itk::Image<unsigned char, DIMENSION> MaskType;
  typedef itk::Point<float, DIMENSION> PointType;
  typedef typename ImageType::OffsetValueType OffsetValueType;

  struct Centroid {
    PointType clInitialPoint, clCurrentPoint;
    float fWeight;

    Centroid() {
      clInitialPoint.Fill(0.0f);
      clCurrentPoint.Fill(0.0f);
      fWeight = 1.0f;
    }

    bool operator<(const Centroid &stOther) const {
      return fWeight < stOther.fWeight;
    }

    bool operator>(const Centroid &stOther) const {
      return fWeight > stOther.fWeight;
    }
  };

  void SetWindow(const PointType &clWindow);

  void SetWindow(float fValue) {
    PointType clWindow;
    clWindow.Fill(fValue);
    SetWindow(clWindow);
  }

  void SetImage(typename ImageType::Pointer p_clImage) {
    Reset();
    m_p_clImage = p_clImage;
  }

  void SetForegroundMask(typename MaskType::Pointer p_clMask) {
    m_p_clMask = p_clMask;
  }

  void Reset() {
    m_p_clImage = NULL;
    m_p_clMask = NULL;
    m_vCentroids.clear();
    m_vNeighbors.clear();
  }

  void AddCentroid(const Centroid &stCentroid) {
    m_vCentroids.push_back(stCentroid);
  }

  void ThresholdCentroids(const PixelType &threshold);

  const std::vector<Centroid> & GetCentroids() const {
    return m_vCentroids;
  }

  void Run();
  void Cluster(unsigned int uiNumCentroids = std::numeric_limits<unsigned int>::max());

private:
  typename ImageType::Pointer m_p_clImage;
  typename MaskType::Pointer m_p_clMask;
  std::vector<Centroid> m_vCentroids;
  std::vector<OffsetValueType> m_vNeighbors;

  bool IsInForeground(const typename MaskType::IndexType &clIndex) const {
    return !m_p_clMask ? true : (m_p_clMask->GetPixel(clIndex) != 0);
  }

  float Update(PointType &clPoint);
};

} // end namespace nih

#include "MeanShift.hpp"

#endif // !MEANSHIFT_H

