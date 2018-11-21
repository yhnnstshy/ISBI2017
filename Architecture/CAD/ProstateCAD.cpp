#include <algorithm>
#include <fstream>
#include <math.h>
#include <limits>
#include <unordered_set>
#include "Common.h"
#include "ProstateCAD.h"
#include "VoiFileIO.h"
#include "VoiToMask.h"
#include <cmath>

namespace nih {

bool ProstateCAD::LoadTrainingList(std::vector<std::string> &vPatientIDs) {
  vPatientIDs.clear();

  const std::string strListFile = GetTrainingListFileName();

  std::ifstream patientStream(strListFile.c_str());
  if (!patientStream) {
    std::cerr << "Error: Could not load patient list file '" << strListFile << "'." << std::endl;
    return false;
  }

  // Patients will be MRN.AccessionNumber
  std::string strLine;
  while (std::getline(patientStream, strLine)) {
    size_t p = strLine.find('\r');

    if (p != std::string::npos)
      strLine.erase(p);

    if (strLine.empty())
      continue;

    vPatientIDs.push_back(strLine);
  }

  return vPatientIDs.size() > 0;
}

bool ProstateCAD::LoadBiopsyData(const char *p_cFileName, std::vector<ProstateCAD::BiopsyPointType> &vBiopsyPoints, const itk::Image<short, 3>::SpacingType &clSpacing) {

  std::ifstream BiopsyFile(p_cFileName);
  //a variable to hold each line of the .csv file. 
  std::string line;
  //variable used to enable string stream. 
  //std::stringstream ss_line;
  //variable used to hold single words or numbers
  //separated by '\t'.
  std::string item;

  int n = 0;
  //Assign x,y,z... values to BiopsyPoint. 
  //Work in progress. Would have to use stringstream.
  //At the same time add the BiopsyPoints into a vector. 
  
  if (BiopsyFile.is_open()){
    while ( std::getline(BiopsyFile, line)){
      std::stringstream ss_line;
      BiopsyPointType b_point;
      ss_line.str(line);
      n = 0;
      while (std::getline(ss_line, item, '\t')){

        if (n == 0)
          b_point.name = item;
        if (n==1) 
          b_point.x = (float)atoi(item.c_str())*clSpacing[0];
        if (n==2)
          b_point.y = (float)atoi(item.c_str())*clSpacing[1];
        if (n==3)
          b_point.z = (float)atoi(item.c_str())*clSpacing[2];
        if (n==4)
          b_point.label = atoi(item.c_str());
        if (n==5)
          b_point.gleason = atoi(item.c_str());
        if (n==6)
          b_point.t2wi = atoi(item.c_str());
        if (n==7)
          b_point.adc = atoi(item.c_str());
        if (n==8)
          b_point.b2000 = atoi(item.c_str());
        n += 1;
      }
      vBiopsyPoints.push_back(b_point);
    }
    BiopsyFile.close();
  }
  else{
    std::cout << "Could not open/find Biopsy file " << std::endl;
    return false;
  }

  return true;
}

bool ProstateCAD::GetVolumeDimensions(const std::string &strPatient, unsigned int a_uiDim[3]) {
  const std::string strPattern = strPatient + ".*.dcm";

  std::vector<std::string> vFileNames;

  FindFiles(vFileNames, GetT2WIDataFolder().c_str(), strPattern.c_str(), false);

  if (vFileNames.empty())
    return false;

  itk::Image<short, 2>::Pointer p_clSlice = LoadDicomSlice(vFileNames[0].c_str());

  if (!p_clSlice) {
    std::cerr << "Error: Could not load slice." << std::endl;
    return false;
  }

  a_uiDim[0] = (unsigned int)p_clSlice->GetBufferedRegion().GetSize()[0];
  a_uiDim[1] = (unsigned int)p_clSlice->GetBufferedRegion().GetSize()[1];
  a_uiDim[2] = (unsigned int)vFileNames.size(); // Number of slices

  return true;
}

bool ProstateCAD::LoadVOIs(const std::string &strFolder, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
  vSlices.clear();

  std::unordered_set<unsigned int> sSlices;

  std::vector<std::string> vVOIFiles;
  FindFiles(vVOIFiles, strFolder.c_str(), "*.voi", false);

  if (vVOIFiles.empty())
    return false;

  for (size_t i = 0; i < vVOIFiles.size(); ++i) {
    const std::string &strFileName = vVOIFiles[i];

    VoiFile clVOIFile;

    if (!clVOIFile.LoadFromFile(strFileName.c_str())) {
      std::cerr << "Error: Could not load '" << strFileName << "'." << std::endl;
      return false;
    }

    const std::vector<VoiFile::Slice> &vVOISlices = clVOIFile.GetSlices();

    for (size_t j = 0; j < vVOISlices.size(); ++j) {
      const VoiFile::Slice &stSlice = vVOISlices[j];
      sSlices.insert(stSlice.uiSliceNumber);
    }

    if (!ConvertVoiToMask3D<unsigned char>(p_clMask, clVOIFile, false)) {
      std::cerr << "Error: Could not convert VOI to mask for file '" << strFileName << "'." << std::endl;
      return false;
    }
  }

  double a_dSpacing[3] = { 1.0, 1.0, 1.0 };

  p_clMask->SetSpacing(a_dSpacing);

  vSlices.insert(vSlices.end(), sSlices.begin(), sSlices.end());
  std::sort(vSlices.begin(), vSlices.end());
  
  return true;
}


bool ProstateCAD::LoadProstate(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
  vSlices.clear();

  std::unordered_set<unsigned int> sSlices;

  const std::string strPatientProstateFile = GetProstateFolder() + '/' + strPatient + ".voi";

  VoiFile clVOIFile;

  if (!clVOIFile.LoadFromFile(strPatientProstateFile.c_str())) {
    std::cerr << "Error: Could not load patient prostate file '" << strPatientProstateFile << "'." << std::endl;
    return false;
  }

  const std::vector<VoiFile::Slice> &vVOISlices = clVOIFile.GetSlices();

  for (size_t i = 0; i < vSlices.size(); ++i) {
    const VoiFile::Slice &stSlice = vVOISlices[i];
    sSlices.insert(stSlice.uiSliceNumber);
  }

  if (!ConvertVoiToMask3D<unsigned char>(p_clMask, clVOIFile, false)) {
    std::cerr << "Error: Could not convert VOI to mask for file '" << strPatientProstateFile << "'." << std::endl;
    return false;
  }

  vSlices.insert(vSlices.end(), sSlices.begin(), sSlices.end());
  std::sort(vSlices.begin(), vSlices.end());

  return true;
}

bool ProstateCAD::LoadCentralGland(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
  vSlices.clear();

  std::unordered_set<unsigned int> sSlices;

  const std::string strPatientCentralGlandFile = GetCentralGlandFolder() + "/voi/" + strPatient + ".voi";

  VoiFile clVOIFile;

  if (!clVOIFile.LoadFromFile(strPatientCentralGlandFile.c_str())) {
    std::cerr << "Error: Could not load patient prostate file '" << strPatientCentralGlandFile << "'." << std::endl;
    return false;
  }

  const std::vector<VoiFile::Slice> &vVOISlices = clVOIFile.GetSlices();

  for (size_t i = 0; i < vSlices.size(); ++i) {
    const VoiFile::Slice &stSlice = vVOISlices[i];
    sSlices.insert(stSlice.uiSliceNumber);
  }

  if (!ConvertVoiToMask3D<unsigned char>(p_clMask, clVOIFile, false)) {
    std::cerr << "Error: Could not convert VOI to mask for file '" << strPatientCentralGlandFile << "'." << std::endl;
    return false;
  }

  vSlices.insert(vSlices.end(), sSlices.begin(), sSlices.end());
  std::sort(vSlices.begin(), vSlices.end());

  return true;
}


} // end namespace nih

