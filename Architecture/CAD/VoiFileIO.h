#ifndef VOIFILEIO_H
#define VOIFILEIO_H

#include <iostream>
#include <fstream>
#include "itkPoint.h"

namespace nih {

class VoiFile {
public:
  typedef itk::Point<float, 2> PointType;

  struct Slice {
    unsigned int uiSliceNumber;
    std::vector<std::vector<PointType> > vContours;

    Slice() {
      uiSliceNumber = 0;
    }

    bool operator<(const Slice &clOther) const {
      return uiSliceNumber < clOther.uiSliceNumber;
    }
  };

  VoiFile() {
    Reset();
  }

  void SaveToStream(std::ostream &os) const;
  bool LoadFromStream(std::istream &is);

  bool SaveToFile(const char *p_cFileName) const {
    std::ofstream fileStream(p_cFileName);

    if (!fileStream)
      return false;

    SaveToStream(fileStream);
    return true;
  }

  bool LoadFromFile(const char *p_cFileName) {
    std::ifstream fileStream(p_cFileName);

    if (!fileStream)
      return false;

    return LoadFromStream(fileStream);
  }

  // As if we ever needed these...
  unsigned int GetCurveType() const {
    return m_uiCurveType;
  }

  void SetCurveType(unsigned int uiCurveType) {
    m_uiCurveType = uiCurveType;
  }

  unsigned int GetColor() const {
    return m_uiColor;
  }

  void SetColor(unsigned int uiColor) {
    m_uiColor = uiColor;
  }

  unsigned int GetUniqueID() const {
    return m_uiUniqueID;
  }

  void SetUniqueID(unsigned int uiUniqueID) {
    m_uiUniqueID = uiUniqueID;
  }

  // But definitely need these
  const std::vector<Slice> & GetSlices() const {
    return m_vSlices;
  }

  std::vector<Slice> & GetSlices() {
    return m_vSlices;
  }

  void Reset() {
    m_uiCurveType = 0;
    m_uiColor = 0;
    m_uiUniqueID = 0;

    m_vSlices.clear();
  }

  bool Empty() const {
    return m_vSlices.empty();
  }

private:
  unsigned int m_uiCurveType, m_uiColor, m_uiUniqueID; // RGBA
  std::vector<Slice> m_vSlices;
};

} // end namespace nih

#endif // !VOIFILEIO_H
