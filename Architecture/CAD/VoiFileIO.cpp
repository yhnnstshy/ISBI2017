#include <ctime>
#include <algorithm>
#include <sstream>
#include <string>
#include <limits>
#include "VoiFileIO.h"

namespace nih {

void VoiFile::SaveToStream(std::ostream &os) const {
  unsigned int uiUniqueID = m_uiUniqueID;

  if (uiUniqueID == 0)
    uiUniqueID = time(NULL);

  os << "MIPAV VOI FILE\n";

  // Curve type
  os << m_uiCurveType << "\t\t# curveType of the VOI\n";

  // Colors (RGBA)
  os << ((m_uiColor >> 24) & 0xff) << "\t\t# color of VOI - red component\n";
  os << ((m_uiColor >> 16) & 0xff) << "\t\t# color of VOI - green component\n";
  os << ((m_uiColor >> 8) & 0xff) << "\t\t# color of VOI - blue component\n";
  os << (m_uiColor & 0xff) << "\t\t# color of VOI - alpha component\n";

  // Number of slices for the VOI
  os << m_vSlices.size() << "\t\t# number of slices for the VOI\n";

  for (size_t i = 0; i < m_vSlices.size(); ++i) {
    const Slice &stSlice = m_vSlices[i];

    // The slice number of the slice
    os << stSlice.uiSliceNumber << "\t\t# slice number\n";

    const std::vector<std::vector<PointType> > &vContours = stSlice.vContours;

    os << vContours.size() << "\t\t# number of contours in slice\n";

    for (size_t j = 0; j < vContours.size(); ++j) {
      const std::vector<PointType> &vContour = vContours[j];

      os << vContour.size() << "\t\t# number of pts in contour\n";

      for (size_t k = 0; k < vContour.size(); ++k) {
        const PointType &clPoint = vContour[k];
        os << clPoint[0] << ' ' << clPoint[1] << '\n';
      }
    }
  }

  os << uiUniqueID << "\t\t# unique ID of the VOI\n";
}

bool VoiFile::LoadFromStream(std::istream &is) {
  Reset();

  is >> std::ws;

  std::string strLine;
  if (!std::getline(is, strLine))
    return false;

  size_t p = strLine.find('\r');
  if (p != std::string::npos)
    strLine.erase(p);

  if (strLine != "MIPAV VOI FILE")
    return false;


  is >> m_uiCurveType;

  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  {
    unsigned int uiComponent = 0;

    if (!(is >> uiComponent) || uiComponent > 255)
      return false;

    // Red
    m_uiColor = (uiComponent << 24);

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (!(is >> uiComponent) || uiComponent > 255)
      return false;

    // Green
    m_uiColor |= (uiComponent << 16);

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (!(is >> uiComponent) ||  uiComponent > 255)
      return false;

    // Blue
    m_uiColor |= (uiComponent << 8);

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (!(is >> uiComponent) || uiComponent > 255)
      return false;

    // Alpha
    m_uiColor |= uiComponent;

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  unsigned int uiNumSlices = 0;

  if (!(is >> uiNumSlices) || uiNumSlices == 0)
    return false;

  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  m_vSlices.resize(uiNumSlices);

  for (size_t i = 0; i < m_vSlices.size(); ++i) {
    Slice &stSlice = m_vSlices[i];

    if (!(is >> stSlice.uiSliceNumber))
      return false;

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    unsigned int uiNumContours = 0;
    if (!(is >> uiNumContours) || uiNumContours == 0)
      return false;

    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::vector<PointType> > &vContours = stSlice.vContours;

    vContours.resize(uiNumContours);

    for (size_t j = 0; j < vContours.size(); ++j) {
      unsigned int uiNumPoints = 0;

      if (!(is >> uiNumPoints) || uiNumPoints == 0)
        return false;

      is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      std::vector<PointType> &vContour = vContours[j];

      vContour.resize(uiNumPoints);

      for (size_t k = 0; k < vContour.size(); ++k) {
        PointType &clPoint = vContour[k];

        if (!(is >> clPoint[0] >> clPoint[1]))
          return false;

        // In case of comments
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }
  }

  m_uiUniqueID = 0;
  if (!(is >> m_uiUniqueID))
    return false;

  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  return true;
}

} // end namespace nih
