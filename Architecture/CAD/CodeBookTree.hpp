#include <cstdint>
#include <string>
#include <queue>
#include "CodeBookTree.h"

template<typename OutputType, typename InputType>
bool SafeAssign(OutputType &y, const InputType &x) {
  if (x < std::numeric_limits<OutputType>::min() || x > std::numeric_limits<OutputType>::max())
    return false;
  y = OutputType(x);
  return true;
}

// CodeBookTreeNode
template<typename InputPixelType, typename OutputPixelType>  
void CodeBookTreeNode<InputPixelType, OutputPixelType>::SaveToStream(std::ostream &os) const {
  os << dThreshold << '\n';
  os << dMinValue << ' ' << dMaxValue << '\n';
  os << (int64_t)code << '\n';
  os << a_uiChildren[0] << ' ' << a_uiChildren[1] << '\n';
}

template<typename InputPixelType, typename OutputPixelType>  
bool CodeBookTreeNode<InputPixelType, OutputPixelType>::LoadFromStream(std::istream &is) {
  if (!(is >> dThreshold))
    return false;

  if (!(is >> dMinValue >> dMaxValue))
    return false;

  if (dMinValue > dMaxValue)
    return false;

  int64_t i64Code = 0;

  if (!(is >> i64Code))
    return false;

  if (!SafeAssign(code, i64Code))
    return false;

  if (!(is >> a_uiChildren[0] >> a_uiChildren[1]))
    return false;

  return true;
}

// CodeBookTree
template<typename InputPixelTypeT, typename OutputPixelTypeT>
int CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::GetIndex(const InputPixelType &pixel) const {
  if (m_vNodes.empty())
    return -1;

  unsigned int uiIndex = 0;

  while (!m_vNodes[uiIndex].IsLeaf()) {
    const NodeType &clNode = m_vNodes[uiIndex];
    uiIndex = clNode.a_uiChildren[(double)pixel > clNode.dThreshold];
  }

  return (int)uiIndex; // XXX: Narrowing
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
bool CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::Invert(const OutputPixelType &inPixel, InputPixelType &outPixel) const {
  if (m_vNodes.empty())
    return false;

  double dMinDist = -1.0;
  double dMean = 0.0;

  for (size_t i = 0; i < m_vNodes.size() && dMinDist != 0.0; ++i) {
    const NodeType &clNode = m_vNodes[i];

    if (!clNode.IsLeaf())
      continue;

    const double dDist = std::abs((double)inPixel - (double)clNode.code);
    if (dMinDist < 0.0 || dDist < dMinDist) {
      dMinDist = dDist;
      dMean = 0.5*(clNode.dMinValue + clNode.dMaxValue);
    }
  }

  if (dMinDist < 0.0)
    return false;

  outPixel = (InputPixelType)(dMean + 0.5f);

  return true;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
bool CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::Train(std::vector<InputPixelType> &vPixels) {
  Clear();

  if (vPixels.size() <= GetMinSampleSize())
    return false;

  if (m_uiMinSampleSize == 0)
    m_uiMinSampleSize = 1;

  m_vNodes.resize(1);
  NodeType &clRoot = m_vNodes[0];

  std::sort(vPixels.begin(), vPixels.end());

  clRoot.p_begin = &vPixels[0];
  clRoot.p_end = clRoot.p_begin + vPixels.size();

  clRoot.dMinValue = (double)vPixels.front();
  clRoot.dMaxValue = (double)vPixels.back();

  std::cout << "Initial KL = " << ComputeKL() << std::endl;

  unsigned int uiLeafCount = 1;
  while (uiLeafCount < GetMaxNumberOfLeaves()) {
    double dMaxGain = -1.0;
    size_t maxNodeIndex = 0;

    for (size_t i = 0; i < m_vNodes.size(); ++i) {
      if (!m_vNodes[i].IsLeaf())
        continue;

      const double dGain = ComputeGain(m_vNodes[i]);
      if (dGain > dMaxGain) {
        dMaxGain = dGain;
        maxNodeIndex = i;
      }
    }

    if (dMaxGain < Small()) // Nothing left to do
      break;

    Grow(m_vNodes[maxNodeIndex]);
    ++uiLeafCount;

    std::cout << "KL = " << ComputeKL() << std::endl;
  }

  std::cout << "Final KL = " << ComputeKL() << std::endl;

  OutputPixelType nextCode = std::numeric_limits<OutputPixelType>::min();
  Label(0, nextCode);

#if 0
  std::cout << "\nIdeal count = " << (size_t)(m_vNodes[0].p_end - m_vNodes[0].p_begin)/(size_t)CountLeaves() << '\n' << std::endl;

  for (size_t i = 0; i < m_vNodes.size(); ++i) {
    const NodeType &clNode = m_vNodes[i];
    if (!clNode.IsLeaf())
      continue;

    std::cout << (int)clNode.code << ' ' << (clNode.p_end - clNode.p_begin) << ' ' << (clNode.p_end - clNode.p_begin)*(clNode.dMaxValue - clNode.dMinValue) << std::endl;
  }
#endif

  return true;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
void CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::SaveToStream(std::ostream &os) const {
  os << GetMagicString() << '\n';
  os << m_vNodes.size() << '\n';

  for (size_t i = 0; i < m_vNodes.size(); ++i)
    m_vNodes[i].SaveToStream(os);
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
bool CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::LoadFromStream(std::istream &is) {
  Clear();

  std::string strMagic;

  is >> std::ws;
  if (!std::getline(is, strMagic))
    return false;

  size_t p = strMagic.find('\r');
  if (p != std::string::npos)
    strMagic.erase(p);

  if (strMagic != GetMagicString())
    return false;

  size_t numNodes = 0;
  if (!(is >> numNodes) || numNodes == 0)
    return false;

  m_vNodes.resize(numNodes);
  for (size_t i = 0; i < m_vNodes.size(); ++i) {
    if (!m_vNodes[i].LoadFromStream(is)) {
      Clear();
      return false;
    }
  }

  return true;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
double CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::ComputeKL() const {
  if (m_vNodes.empty())
    return -1.0;

  double dNormalizer = 0.0;

  const size_t totalSampleSize = m_vNodes[0].p_end - m_vNodes[0].p_begin;
  const double dGTProb = 1.0 / (m_vNodes[0].dMaxValue - m_vNodes[0].dMinValue);

  for (size_t i = 0; i < m_vNodes.size(); ++i) {
    const NodeType &clNode = m_vNodes[i];

    if (!clNode.IsLeaf())
      continue;

    const size_t sampleSize = clNode.p_end - clNode.p_begin;
    const double dLength = clNode.dMaxValue - clNode.dMinValue;

    const double dProb = (double)sampleSize / (double)totalSampleSize;

    dNormalizer += dLength * dProb;
  }

  double dKL = 0.0;
  for (size_t i = 0; i < m_vNodes.size(); ++i) {
    const NodeType &clNode = m_vNodes[i];

    if (!clNode.IsLeaf())
      continue;

    const size_t sampleSize = clNode.p_end - clNode.p_begin;
    const double dLength = clNode.dMaxValue - clNode.dMinValue;

    double dProb = (double)sampleSize / (double)totalSampleSize;
    dProb /= dNormalizer;

    dKL += dLength * dGTProb * std::log(dGTProb / dProb);
  }

  return dKL;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
double CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::ComputeGain(NodeType &clNode) const {
  if (clNode.dMaxGain >= 0)
    return clNode.dMaxGain;

  const size_t sampleSize = clNode.p_end - clNode.p_begin;
  if (sampleSize < GetMinSampleSize() || clNode.dMaxValue - clNode.dMinValue <= 1.0)
    return -1.0;

  std::vector<double> vThresholds;
  ComputeThresholds(clNode, vThresholds);

  if (vThresholds.empty())
    return -1.0;

  std::vector<size_t> vBinCount(vThresholds.size(), 0);

  for (const InputPixelType *itr = clNode.p_begin; itr != clNode.p_end; ++itr) {
    const double dValue = (double)(*itr);

    for (size_t i = 0; i < vThresholds.size() && dValue > vThresholds[i]; ++i)
      ++vBinCount[i];
  }

  for (size_t i = 0; i < vThresholds.size(); ++i) {
    if (vBinCount[i] == 0 || vBinCount[i] == sampleSize)
      continue;

    const double dThreshold = vThresholds[i];

    const double dRightProb = (double)vBinCount[i] / (double)sampleSize;
    const double dRightLength = clNode.dMaxValue - dThreshold;

    const double dLeftProb = (double)(sampleSize - vBinCount[i])/(double)sampleSize;
    const double dLeftLength = dThreshold - clNode.dMinValue;

    const double dGain = -dLeftLength * std::log(dLeftProb) - dRightLength * std::log(dRightProb);

    if (dGain > clNode.dMaxGain) {
      clNode.dMaxGain = dGain;
      clNode.dThreshold = dThreshold;
    }
  }

  return clNode.dMaxGain;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
void CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::ComputeThresholds(NodeType &clNode, std::vector<double> &vThresholds) const {
  vThresholds.clear();

  double dPrevValue = (double)(*clNode.p_begin);

  for (const InputPixelType *itr = clNode.p_begin+1; itr != clNode.p_end; ++itr) {
    const double dValue = (double)(*itr);

    if (dValue - dPrevValue > Small()) {
      const double dThreshold = 0.5*(dValue + dPrevValue);
      vThresholds.push_back(dThreshold);
    }

    dPrevValue = dValue;
  }
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
bool CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::Grow(NodeType &clNode) {
  const InputPixelType *splitItr = NULL;

  for (splitItr = clNode.p_begin; splitItr != clNode.p_end && (double)(*splitItr) <= clNode.dThreshold; ++splitItr);

  if (splitItr == clNode.p_begin || splitItr == clNode.p_end)
    return false;

  clNode.a_uiChildren[0] = (unsigned int)m_vNodes.size();
  clNode.a_uiChildren[1] = clNode.a_uiChildren[0] + 1;

  NodeType a_clChildren[2];

  a_clChildren[0].p_begin = clNode.p_begin;
  a_clChildren[0].p_end = splitItr;

  a_clChildren[0].dMinValue = clNode.dMinValue;
  a_clChildren[0].dMaxValue = clNode.dThreshold;

  a_clChildren[1].p_begin = splitItr;
  a_clChildren[1].p_end = clNode.p_end;

  a_clChildren[1].dMinValue = clNode.dThreshold;
  a_clChildren[1].dMaxValue = clNode.dMaxValue;

  m_vNodes.push_back(a_clChildren[0]);
  m_vNodes.push_back(a_clChildren[1]);

  return true;
}

template<typename InputPixelTypeT, typename OutputPixelTypeT>
void CodeBookTree<InputPixelTypeT,OutputPixelTypeT,true>::Label(unsigned int uiIndex, OutputPixelType &nextCode) {
  NodeType &clNode = m_vNodes[uiIndex];

  // In-order visit to keep codes consistently ordered with InputPixelType

  if (clNode.IsLeaf()) {
    clNode.code = nextCode++;
  }
  else {
    Label(clNode.a_uiChildren[0], nextCode);
    // Nothing to do at this node
    Label(clNode.a_uiChildren[1], nextCode);
  }
}
