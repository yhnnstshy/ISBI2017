#include <cmath>
#include <algorithm>
#include <unordered_map>
#include "MeanShift.h"

namespace nih {

namespace MeanShiftUtils {

template<unsigned int Index, unsigned int Dimension>
struct Incrementer {
  typedef itk::Point<float, Dimension> PointType;

  static void Run(PointType &clPoint, const PointType &clBegin, const PointType &clEnd) {
    clPoint[Index] += 1.0f;
    if (Index+1 < Dimension && clPoint[Index] >= clEnd[Index]) {
      clPoint[Index] = clBegin[Index];
      Incrementer<Index+1,Dimension>::Run(clPoint,clBegin,clEnd);
    }
  }
};

template<unsigned int Dimension>
struct Incrementer<Dimension,Dimension> {
  typedef itk::Point<float, Dimension> PointType;

  static void Run(PointType &clPoint, const PointType &clBegin, const PointType &clEnd) { }
};

template<unsigned int Dimension>
void Increment(itk::Point<float, Dimension> &clPoint, const itk::Point<float, Dimension> &clBegin, const itk::Point<float, Dimension> &clEnd) {
  Incrementer<0, Dimension>::Run(clPoint, clBegin, clEnd);
}

} // end namespace MeanShiftUtils

template<typename PixelType, unsigned int Dimension>
void MeanShift<PixelType, Dimension>::SetWindow(const PointType &clWindow) {
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;

  if (!m_p_clImage)
    return;

  const SizeType &clSize = m_p_clImage->GetBufferedRegion().GetSize();

  m_vNeighbors.clear();

  PointType clBegin, clEnd;
  for (int i = 0; i < DIMENSION; ++i) {
    const float fWindowSize = std::min((float)clSize[i], clWindow[i]);

    clBegin[i] = -(int)(fWindowSize/2);
    clEnd[i] = clBegin[i] + fWindowSize;
  }
  
  PointType clCurrent = clBegin;
  while (clCurrent[DIMENSION-1] != clEnd[DIMENSION-1]) {

    //std::cout << clCurrent << std::endl;

    IndexType clIndex;

    float r = 0.0f;
    for (int i = 0; i < DIMENSION; ++i)
      r += std::pow(clCurrent[i]/clWindow[i], 2);

    clIndex.CopyWithCast(clCurrent);

    if (r < 1.0f) {
      OffsetValueType offset = m_p_clImage->ComputeOffset(clIndex);

      m_vNeighbors.push_back(offset);
    }

    MeanShiftUtils::Increment(clCurrent, clBegin, clEnd);
  }

  // Add the last one
  {
    IndexType clIndex;

    float r = 0.0f;
    for (int i = 0; i < DIMENSION; ++i)
      r += std::pow(clCurrent[i]/clWindow[i], 2);

    if (r < 1.0f) {
      OffsetValueType offset = m_p_clImage->ComputeOffset(clIndex);

      m_vNeighbors.push_back(offset);
    }
  }
}

template<typename PixelType, unsigned int Dimension>
void MeanShift<PixelType, Dimension>::ThresholdCentroids(const PixelType &threshold) {
  typedef typename ImageType::IndexType IndexType;

  if (!m_p_clImage) {
    return;
  }

  const PixelType * const p_buffer = m_p_clImage->GetPixelContainer()->GetImportPointer();
  const size_t bufferSize = m_p_clImage->GetPixelContainer()->Size();

  if (p_buffer == NULL)
    return;

  Centroid stCentroid;
  stCentroid.fWeight = 1.0f;

  for (size_t i = 0; i < bufferSize; ++i) {
    if (p_buffer[i] > threshold) {
      IndexType clIndex = m_p_clImage->ComputeIndex(i);

      if (!IsInForeground(clIndex))
        continue;

      for (int j = 0; j < DIMENSION; ++j)
        stCentroid.clInitialPoint[j] = (float)clIndex[j];

      m_vCentroids.push_back(stCentroid);
    }
  }

}

template<typename PixelType, unsigned int Dimension>
void MeanShift<PixelType, Dimension>::Run() {
  const int iMaxSteps = 1000000;
  //const float fTol = 1e-6f;
  const float fTol = 0.0f;

  //std::cout << "Number of centroids: " << m_vCentroids.size() << std::endl;

  for (size_t i = 0; i < m_vCentroids.size(); ++i)
    m_vCentroids[i].clCurrentPoint = m_vCentroids[i].clInitialPoint;

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; (size_t)i < m_vCentroids.size(); ++i) {
    float fUpdateNorm = fTol + 1.0f;

    for (int s = 0; s < iMaxSteps && fUpdateNorm > fTol; ++s) {
      PointType clBefore = m_vCentroids[i].clCurrentPoint;

      fUpdateNorm = Update(m_vCentroids[i].clCurrentPoint);

//#pragma omp critical
      //std::cout << clBefore << " -> " << m_vCentroids[i].clCurrentPoint << ": norm = " << fUpdateNorm << std::endl;
    }
  }
}

template<typename PixelType, unsigned int Dimension>
void MeanShift<PixelType, Dimension>::Cluster(unsigned int uiNumCentroids) {
  typedef typename ImageType::IndexType IndexType;
  typedef std::unordered_map<OffsetValueType, Centroid> MapType;

  if (!m_p_clImage || m_vCentroids.empty())
    return;

#if 1
  typedef itk::Image<unsigned char, DIMENSION> MaskType;
  typename MaskType::Pointer p_clMask = MaskType::New();

  p_clMask->SetRegions(m_p_clImage->GetBufferedRegion());
  p_clMask->SetSpacing(m_p_clImage->GetSpacing());
  p_clMask->SetOrigin(m_p_clImage->GetOrigin());

  p_clMask->Allocate(true);

  for (size_t i = 0; i < m_vCentroids.size(); ++i) {
    const float fWeight = m_vCentroids[i].fWeight;
    IndexType clIndex;
    clIndex.CopyWithCast(m_vCentroids[i].clCurrentPoint);

    unsigned int uiVoxel = p_clMask->GetPixel(clIndex);
    uiVoxel += (int)fWeight;

    uiVoxel = std::min(uiVoxel, (unsigned int)255);

    p_clMask->SetPixel(clIndex, (unsigned char)uiVoxel);
  }

  nih::SaveVolume<unsigned char>(p_clMask, "MeanShiftMask.mhd");
#endif

  MapType mNewCentroids;
  typename MapType::iterator itr;

  for (size_t i = 0; i < m_vCentroids.size(); ++i) {
    const float fWeight = m_vCentroids[i].fWeight;

    IndexType clIndex;
    clIndex.CopyWithCast(m_vCentroids[i].clCurrentPoint);

    const OffsetValueType index = m_p_clImage->ComputeOffset(clIndex);
    itr = mNewCentroids.find(index);

    if (itr == mNewCentroids.end()) {
      Centroid &stNewCentroid = mNewCentroids[index];
      stNewCentroid.fWeight = fWeight;
      stNewCentroid.clCurrentPoint = m_vCentroids[i].clCurrentPoint;
    }
    else {
      Centroid &stNewCentroid = itr->second;
      stNewCentroid.fWeight += fWeight;
    }
  }

  m_vCentroids.clear();
  for (itr = mNewCentroids.begin(); itr != mNewCentroids.end(); ++itr)
    m_vCentroids.push_back(itr->second);

  std::sort(m_vCentroids.begin(), m_vCentroids.end(), std::greater<Centroid>());

  size_t newSize = std::min((size_t)uiNumCentroids, m_vCentroids.size());

  m_vCentroids.resize(newSize);
}

template<typename PixelType, unsigned int Dimension>
float MeanShift<PixelType, Dimension>::Update(PointType &clPoint) {
  typedef typename ImageType::IndexType IndexType;

  IndexType clIndex;

  clIndex.CopyWithCast(clPoint);

  const OffsetValueType index = m_p_clImage->ComputeOffset(clIndex);

  if (!m_p_clImage->GetBufferedRegion().IsInside(clIndex))
    return -1.0f;

  PointType clNewPoint;
  clNewPoint.Fill(0.0f);

  float fWeightSum = 0.0f;

  for (size_t i = 0; i < m_vNeighbors.size(); ++i) {
    const OffsetValueType neighborIndex = index + m_vNeighbors[i];

    clIndex = m_p_clImage->ComputeIndex(neighborIndex);

    if (m_p_clImage->GetBufferedRegion().IsInside(clIndex) && IsInForeground(clIndex)) {
      const float fWeight = (float)(m_p_clImage->GetPixel(clIndex));
      fWeightSum += fWeight;

      for (int j = 0; j < DIMENSION; ++j)
        clNewPoint[j] += (float)(fWeight * clIndex[j]);
    }
  }

  float fDist2 = 0.0f;

  if (fWeightSum > 0.0f) {
    for (int j = 0; j < DIMENSION; ++j)
      clNewPoint[j] /= fWeightSum;

    fDist2 = clPoint.SquaredEuclideanDistanceTo(clNewPoint);

    clPoint = clNewPoint;
  }


  return fDist2;
}

} // end namespace nih
