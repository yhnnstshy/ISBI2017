#ifndef PARAMETERCONTAINER_H
#define PARAMETERCONTAINER_H

#include <string>
#include <sstream>

namespace nih {
	
class ParameterContainer {
public:
  virtual ~ParameterContainer() { }

  template<typename Tp>
  Tp GetValue(const char *pKey, const Tp &defaultValue) const;

  template<typename Tp>
  void SetValue(const char *pKey, const Tp &value);

  bool HasKey(const char *pKey) const {
    std::string strValue;
    return GetStringValue(pKey, strValue);
  }

protected:
  virtual bool GetStringValue(const char *pKey, std::string &strValue) const = 0;
  virtual void SetStringValue(const char *pKey, const std::string &strValue) = 0;
};

template<typename Tp>
Tp ParameterContainer::GetValue(const char *pKey, const Tp &defaultValue) const {
  std::string strValue;
  
  if (!GetStringValue(pKey, strValue))
    return defaultValue;

  std::stringstream convertStream;
  convertStream.str(strValue);

  Tp value(defaultValue);

  return !(convertStream >> value) ? defaultValue : value;
}

template<>
inline std::string ParameterContainer::GetValue<std::string>(const char *pKey, const std::string &defaultValue) const {
  std::string strValue;

  if (!GetStringValue(pKey, strValue))
    return defaultValue;

  return strValue;
}

template<>
inline bool ParameterContainer::GetValue<bool>(const char *pKey, const bool &bDefaultValue) const {
  std::string strValue;

  if (!GetStringValue(pKey, strValue))
    return bDefaultValue;

  std::stringstream convertStream;
  convertStream.str(strValue);

  bool bValue(bDefaultValue);

  // Try it the usual way first
  if (!(convertStream >> bValue)) {
    convertStream.clear();
    convertStream.str(strValue);
  }
  else // Successfully read
    return bValue;

  return !(convertStream >> std::boolalpha >> bValue) ? bDefaultValue : bValue;
}

template<typename Tp>
void ParameterContainer::SetValue(const char *pKey, const Tp &value) {
  std::stringstream convertStream;

  convertStream << value;

  SetStringValue(pKey, convertStream.str());
}

template<>
inline void ParameterContainer::SetValue<std::string>(const char *pKey, const std::string &strValue) {
  SetStringValue(pKey, strValue);
}

template<>
inline void ParameterContainer::SetValue<bool>(const char *pKey, const bool &bValue) {
  std::stringstream convertStream;

  convertStream << std::boolalpha << bValue;

  SetStringValue(pKey, convertStream.str());
}

} // end namespace nih

#endif // !PARAMETERCONTAINER_H

