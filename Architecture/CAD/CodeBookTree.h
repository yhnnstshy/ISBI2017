#ifndef CODEBOOKTREE_H
#define CODEBOOKTREE_H

#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <type_traits>

template<typename InputPixelType, typename OutputPixelType>
class CodeBookTreeNode {
public:
  typedef CodeBookTreeNode<InputPixelType, OutputPixelType> Self;

  OutputPixelType code;
  double dThreshold, dMinValue, dMaxValue;
  unsigned int a_uiChildren[2];

  const InputPixelType *p_begin, *p_end;
  double dMaxGain;

  CodeBookTreeNode() {
    code = OutputPixelType();
    dThreshold = dMinValue = dMaxValue = 0.0;
    a_uiChildren[0] = a_uiChildren[1] = 0;
    p_begin = p_end = NULL;
    dMaxGain = -1.0;
  }

  bool IsLeaf() const {
    return a_uiChildren[0] == 0;
  }

  void SaveToStream(std::ostream &os) const;
  bool LoadFromStream(std::istream &is);
};

template<typename InputPixelTypeT, typename OutputPixelTypeT, bool Integral = std::is_integral<OutputPixelTypeT>::value>
class CodeBookTree { };

template<typename InputPixelTypeT, typename OutputPixelTypeT>
class CodeBookTree<InputPixelTypeT, OutputPixelTypeT, true> {
public:
  typedef InputPixelTypeT InputPixelType;
  typedef OutputPixelTypeT OutputPixelType;
  typedef CodeBookTreeNode<InputPixelType, OutputPixelType> NodeType;

  static unsigned int GetMaxNumberOfLeaves() {
    return (unsigned int)std::numeric_limits<OutputPixelType>::max() - (unsigned int)std::numeric_limits<OutputPixelType>::min() + 1;
  }

  static double Small() {
    return 1e-7;
  }

  static const char * GetMagicString() {
    return "CODEBOOKTREE";
  }

  CodeBookTree() {
    m_uiMinSampleSize = 100;
  }

  unsigned int GetMinSampleSize() const {
    return m_uiMinSampleSize;
  }

  void SetMinSampleSize(unsigned int uiMinSampleSize) {
    m_uiMinSampleSize = uiMinSampleSize;
  }

  void Clear() {
    m_vNodes.clear();
  }

  const std::vector<NodeType> & GetNodes() const {
    return m_vNodes;
  }

  bool Empty() const {
    return m_vNodes.empty();
  }

  unsigned int CountLeaves() const {
    unsigned int uiCount = 0;
    for (size_t i = 0; i < m_vNodes.size(); ++i) {
      if (m_vNodes[i].IsLeaf()) {
        ++uiCount;
      }
    }

    return uiCount;
  }

  InputPixelType MinPixel() const {
    return m_vNodes.empty() ? InputPixelType() : InputPixelType(m_vNodes[0].dMinValue);
  }

  InputPixelType MaxPixel() const {
    return m_vNodes.empty() ? InputPixelType() : InputPixelType(m_vNodes[0].dMaxValue);
  }

  int GetIndex(const InputPixelType &pixel) const;

  bool Convert(const InputPixelType &inPixel, OutputPixelType &outPixel) const {
    const int iIndex = GetIndex(inPixel);
    if (iIndex < 0)
      return false;

    outPixel = m_vNodes[iIndex].code;

    return true;
  }

  bool Invert(const OutputPixelType &inPixel, InputPixelType &outPixel) const;

  bool Train(std::vector<InputPixelType> &vPixels);

  void SaveToStream(std::ostream &os) const;
  bool LoadFromStream(std::istream &is);

  bool SaveToFile(const char *p_cFileName) const {
    std::ofstream outStream(p_cFileName);
    if (!outStream)
      return false;

    SaveToStream(outStream);
    return true;
  }

  bool LoadFromFile(const char *p_cFileName) {
    std::ifstream inStream(p_cFileName);
    if (!inStream)
      return false;

    return LoadFromStream(inStream);
  }

private:
  unsigned int m_uiMinSampleSize;
  std::vector<NodeType> m_vNodes;

  double ComputeKL() const;

  double ComputeGain(NodeType &clNode) const;
  void ComputeThresholds(NodeType &clNode, std::vector<double> &vThresholds) const;
  bool Grow(NodeType &clNode);
  void Label(unsigned int uiIndex, OutputPixelType &nextCode);
};

#include "CodeBookTree.hpp"

#endif // !CODEBOOKTREE_H
