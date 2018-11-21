#ifndef PARAMETERFILE_H
#define PARAMETERFILE_H

#include "ParameterContainer.h"
#include "IniFileIO.h"

namespace nih {

class ParameterFile : public ParameterContainer {
public:
  ParameterFile(const char *pFileName, const char *pINIGroup) {
    m_strFileName = pFileName;
    m_strINIGroup = pINIGroup;
  }

  virtual ~ParameterFile() { }

  virtual void SetGroup(const char *pINIGroup) {
    m_strINIGroup = pINIGroup;
  }

protected:
  virtual bool GetStringValue(const char *pKey, std::string &strValue) const {
    strValue = IniFileIO::ReadProfileString(m_strFileName, m_strINIGroup, pKey, "DEFAULT_VALUE");
    return strValue != "DEFAULT_VALUE";
  }

  virtual void SetStringValue(const char *pKey, const std::string &strValue) {
    IniFileIO::WriteProfileString(m_strFileName, m_strINIGroup, pKey, strValue);
  }

private:
  std::string m_strFileName;
  std::string m_strINIGroup;
};

} // end namespace nih

#endif // !PARAMETERFILE_H

