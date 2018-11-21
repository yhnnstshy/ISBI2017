#ifndef PROSTATECAD_H
#define PROSTATECAD_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "itkImage.h"

namespace nih {

class ProstateCAD {
public:
  enum { DIMENSION = 2 };
  //BiopysPoint struct start. 
  struct BiopsyPoint{
    std::string name;
    float x, y, z;
    int label;
    int gleason;
    int t2wi;
    int adc;
    int b2000;

    BiopsyPoint() {
      x = y = z = 0.0f;
      label = gleason = t2wi = adc = b2000 = 0;
    }
  };
  typedef BiopsyPoint BiopsyPointType;

  ProstateCAD() {
    SetModelFolder(".");
    SetDataFolder(".");
  }

  virtual ~ProstateCAD() { }

  void SetModelFolder(const std::string &strModelFolder) {
    m_strModelFolder = strModelFolder;
  }

  void SetDataFolder(const std::string &strDataFolder) {
    m_strDataFolder = strDataFolder;
  }

  std::string GetModelFolder() const {
    return m_strModelFolder;
  }

  std::string GetDataFolder() const {
    return m_strDataFolder;
  }

  // T2WI volume data folder
  std::string GetT2WIDataFolder() const {
    return GetDataFolder() + "/T2WI";
  }

  // sT2Map volume data folder
  std::string GetsT2MapDataFolder() const {
    return GetDataFolder() + "/sT2Map";
  }


  // ADC volume data folder
  std::string GetADCDataFolder() const {
    return GetDataFolder() + "/ADC";
  }

  // B2000 volume data folder
  std::string GetB2000DataFolder() const {
    return GetDataFolder() + "/B2000";
  }

  // In each folder are a set of VOI files (with any name) with the tumors. Each file is one or more tumors (usually just one).
  std::string GetTumorFolder() const {
    return GetDataFolder() + "/Tumors";
  }

  // Similar to GetTumorFolder()
  std::string GetBPHFolder() const {
    return GetDataFolder() + "/BPHs";
  }

  // VOI prostate segmentation folder
  // Folder contains files MRN.AccessionNumber.voi
  std::string GetProstateFolder() const {
    return GetDataFolder() + "/Prostates";
  }
  std::string GetBiopsyFile() const {
    return GetDataFolder() + "/BiopsyData.csv";
  }

  // VOI central gland segmentation folder
  // Folder contains files MRN.AccessionNumber.voi
  std::string GetCentralGlandFolder() const {
    return GetDataFolder() + "/CentralGlands";
  }


  std::string GetTrainingListFileName() const {
    return GetModelFolder() + "/PatientList.txt";
  }


  bool LoadTrainingList(std::vector<std::string> &vPatientIDs);


  bool GetVolumeDimensions(const std::string &strPatient, unsigned int a_uiDim[3]);

  bool LoadVOIs(const std::string &strFolder, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices);

  // This function loads all tumors for a given patient (MRN.AccessionNumber) and returns a 3D binary mask for the tumors as well as the slices where tumors exist
  bool LoadTumors(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
    return LoadVOIs(GetTumorFolder() + '/' + strPatient, p_clMask, vSlices);
  }

  // This function loads all BPHs for a given patient (MRN.AccessionNumber) and returns a 3D binary mask for the tumors as well as the slices where tumors exist
  bool LoadBPHs(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
    return LoadVOIs(GetBPHFolder() + '/' + strPatient, p_clMask, vSlices);
  }

  // This function loads the prostate segmentation for a given patient (MRN.AccessionNumber) and returns a 3D binary mask for the prostate as well as the slices where the prostate exists.
  bool LoadProstate(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices);
  bool LoadCentralGland(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices);

  static bool LoadBiopsyData(const char *p_cFileName, std::vector<ProstateCAD::BiopsyPointType> &vBiopsyPoints, const itk::Image<short, 3>::SpacingType &clSpacing);
  bool LoadBiopsyData(std::vector<ProstateCAD::BiopsyPointType> &vBiopsyPoints, const itk::Image<short, 3>::SpacingType &clSpacing) {
    return LoadBiopsyData(GetBiopsyFile().c_str(), vBiopsyPoints, clSpacing);
  }

  virtual void SetForeground(itk::Image<unsigned char, DIMENSION>::Pointer p_clForeground) {
    m_p_clForeground = p_clForeground;
  }


protected:
  float m_frmin, m_frmid, m_frmax;

  std::string m_strDataFolder, m_strModelFolder;

  itk::Image<unsigned char, DIMENSION>::Pointer m_p_clForeground;

  virtual bool Foreground(int x, int y) const {
    if (!m_p_clForeground)
      return true;

    const itk::Index<2> clIndex = { x, y };
    return m_p_clForeground->GetPixel(clIndex) != 0;
  }

};

} // end namespace nih

#endif // !PROSTATECAD_H
