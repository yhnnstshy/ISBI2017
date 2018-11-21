#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <map>
#include "ParameterContainer.h"

namespace nih {

class ParameterMap : public ParameterContainer {
public:
  virtual ~ParameterMap() { }

protected:
  virtual bool GetStringValue(const char *pKey, std::string &strValue) const {
    std::map<std::string, std::string>::const_iterator itr = m_mKeyValueMap.find(pKey);

    if (itr == m_mKeyValueMap.end())
      return false;

    strValue = itr->second;

    return true;
  }

  virtual void SetStringValue(const char *pKey, const std::string &strValue) {
    m_mKeyValueMap[pKey] = strValue;
  }

private:
  std::map<std::string, std::string> m_mKeyValueMap;
};

} // end namespace nih

#endif // !PARAMETERMAP_H

