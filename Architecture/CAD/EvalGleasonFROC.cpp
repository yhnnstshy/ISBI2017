#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "ProstateCAD.h"
#include "Common.h"
#include "bsdgetopt.h"

class DetectionHistogram {
public:
  enum { NUM_THRESHOLDS = 1000 };

  DetectionHistogram() {
    Clear();
  }  

  void Clear() {
    std::fill(m_a_uiBinCount, m_a_uiBinCount + GetNumBins(), (unsigned int)0);
    m_uiPatientCount = m_uiTotalCount = 0;
  }

  unsigned int GetNumBins() const {
    return (unsigned int)NUM_THRESHOLDS;
  }
  
  unsigned int GetTotalCount() const {
    return m_uiTotalCount;
  }

  unsigned int GetPatientCount() const {
    return m_uiPatientCount;
  }

  float GetThresholdByIndex(unsigned int uiBin) const {
    const float fDelta = (GetMax() - GetMin()) / (float)GetNumBins();
    return GetMin() + fDelta * (uiBin + 0.5f);
  }

  unsigned int GetCountByIndex(unsigned int uiBin) const {
    return uiBin < GetNumBins() ? m_a_uiBinCount[uiBin] : 0;
  }

  float GetRateByIndex(unsigned int uiBin) const {
    return m_uiTotalCount > 0 ? GetCountByIndex(uiBin) / (float)m_uiTotalCount : 0.0f;
  }

  float GetAverageByIndex(unsigned int uiBin) const {
    const unsigned int uiPatientCount = GetPatientCount();
    return uiPatientCount > 0 ? GetCountByIndex(uiBin) / (float)uiPatientCount : 0.0f;
  }

  unsigned int GetCountByThreshold(float fThreshold) const {
    if (fThreshold < GetMin() || fThreshold > GetMax())
      return 0;

    unsigned int uiBin = (unsigned int)( GetNumBins() * (fThreshold - GetMin())/(GetMax() - GetMin()) + 0.5f );
    uiBin = std::min(GetNumBins()-1, uiBin);

    return m_a_uiBinCount[uiBin];
  }

  float GetMin() const {
    return 0.0;
  }

  float GetMax() const {
    return 1.0;
  }

  void Update(float fValue) {
    const float fDelta = (GetMax() - GetMin()) / (float)GetNumBins();

    ++m_uiTotalCount;

    for (unsigned int i = 0; i < GetNumBins() && fValue > GetMin() + fDelta*(i + 0.5f); ++i)
      ++m_a_uiBinCount[i];
  }

  void NextPatient() {
    ++m_uiPatientCount;
  }

  void SetPatientCount(unsigned int uiCount) {
    m_uiPatientCount = uiCount;
  }

  void Merge(const DetectionHistogram &clOther) {
    for (unsigned int i = 0; i < GetNumBins(); ++i)
      m_a_uiBinCount[i] += clOther.m_a_uiBinCount[i];

    m_uiPatientCount += clOther.m_uiPatientCount;
    m_uiTotalCount += clOther.m_uiTotalCount;
  }

private:
  unsigned int m_a_uiBinCount[NUM_THRESHOLDS];
  unsigned int m_uiPatientCount, m_uiTotalCount;
};

class BindFROC {
public:
  const DetectionHistogram &clTruePositives, &clFalsePositives;

  BindFROC(const DetectionHistogram &clTruePositives_, const DetectionHistogram &clFalsePositives_) 
  : clTruePositives(clTruePositives_), clFalsePositives(clFalsePositives_) { }

  void Print(std::ostream &os) const {
    os << "# Threshold\tTrue Positive Rate\tAverage False Positives\n";
    for (unsigned int i = 0; i < clTruePositives.GetNumBins(); ++i)
      os << clTruePositives.GetThresholdByIndex(i) << '\t' << clTruePositives.GetRateByIndex(i) << '\t' << clFalsePositives.GetAverageByIndex(i) << '\n';
  }
};

class BindROC {
public:
  const DetectionHistogram &clTruePositives, &clFalsePositives;

  BindROC(const DetectionHistogram &clTruePositives_, const DetectionHistogram &clFalsePositives_) 
  : clTruePositives(clTruePositives_), clFalsePositives(clFalsePositives_) { }

  void Print(std::ostream &os) const {
    os << "# Threshold\tTrue Positive Rate\tFalse Positive Rate\n";
    for (unsigned int i = 0; i < clTruePositives.GetNumBins(); ++i)
      os << clTruePositives.GetThresholdByIndex(i) << '\t' << clTruePositives.GetRateByIndex(i) << '\t' << clFalsePositives.GetRateByIndex(i) << '\n';
  }
};

std::ostream & operator<<(std::ostream &os, const BindFROC &clBind) {
  clBind.Print(os);
  return os;
}

std::ostream & operator<<(std::ostream &os, const BindROC &clBind) {
  clBind.Print(os);
  return os;
}

class DetectionGraph {
public:
  typedef itk::Point<float, 3> PointType;

  DetectionGraph(const std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints, const std::vector<std::pair<PointType, float> > &vDetectionPairs) {
    m_vBiopsyPoints = vBiopsyPoints;
    m_vDetectionPairs = vDetectionPairs;

    ComputeGraph();
  }

  float DistanceThreshold() const { return 10.0f; }

  float Optimize(std::vector<int> &vPairs) const {
    vPairs.resize(m_vBiopsyPoints.size());
    std::fill(vPairs.begin(), vPairs.end(), -1);

    std::vector<int> vClaimed(m_vDetectionPairs.size(), 0);

    return OptimizeHelper(0, 0.0f, vPairs, vClaimed);
  }

private:
  std::vector<nih::ProstateCAD::BiopsyPointType> m_vBiopsyPoints;
  std::vector<std::pair<PointType, float> > m_vDetectionPairs;

  std::vector<std::vector<int> > m_vGraph; // m_vGraph[biopsyPoint] --> List of neighbors

  float Distance(float x1, float y1, float z1, float x2, float y2, float z2) const {
    return std::sqrt(std::pow(x1-x2, 2) + std::pow(y1-y2, 2) + std::pow(z1-z2, 2));
  }

  float OptimizeHelper(int iBiopsyPoint, float fProbSum, std::vector<int> &vPairs, std::vector<int> &vClaimed) const;

  void ComputeGraph();
};

class DetectionKernel {
public:
  DetectionKernel(unsigned int uiWidth, unsigned int uiHeight) {
    m_uiWidth = uiWidth;
    m_uiHeight = uiHeight;

    const float fSigma = m_uiWidth * m_uiHeight / (-4.0f * std::log(0.01f));
    m_vKernel.resize(m_uiWidth * m_uiHeight);

    for (unsigned int j = 0; j < uiHeight; ++j) {
      const float fY = j - 0.5f*uiHeight;

      for (unsigned int i = 0; i < uiWidth; ++i) {
        const float fX = i - 0.5f*uiWidth;

        m_vKernel[uiWidth*j + i] = std::exp(-(fX*fX + fY*fY)/fSigma);
      }
    }

    const float fSum = std::accumulate(m_vKernel.begin(), m_vKernel.end(), 0.0f);

    if (fSum <= 0)
      return;

    for (size_t i = 0; i < m_vKernel.size(); ++i)
      m_vKernel[i] /= fSum;
  }

  float Evaluate(int x, int y, int z, itk::Image<short, 3>::Pointer p_clProbMap) {
    typedef itk::Image<short, 3> ImageType;
    const ImageType::SizeType &clSize = p_clProbMap->GetBufferedRegion().GetSize();

    const int yBegin = y - (int)m_uiHeight/2;
    const int yEnd = yBegin + m_uiHeight;

    const int xBegin = x - (int)m_uiWidth/2;
    const int xEnd = xBegin + m_uiWidth;

    float fSum = 0.0f;

    for (int j = yBegin; j < yEnd; ++j) {
      if (j < 0 || j >= clSize[1])
        continue;

      for (int i = xBegin; i < xEnd; ++i) {
        if (i < 0 || i >= clSize[0])
          continue;

        const ImageType::IndexType clIndex = { i, j, z };
        const unsigned int uiIndex = m_uiWidth * (j - yBegin) + (i - xBegin);

        fSum += m_vKernel[uiIndex] * p_clProbMap->GetPixel(clIndex);
      }
    }

    return fSum;
  }

private:
  unsigned int m_uiWidth, m_uiHeight;
  std::vector<float> m_vKernel;
};

void Usage(const char *p_cArg0) {
  std::cerr << "Usage: " << p_cArg0 << " [-h] -d dataFolder -l biopsyList -n normalList -p probabilityMapFolder" << std::endl;
  exit(1);
}

class EvalGleasonFROC {
public:
  enum { NUM_THRESHOLDS = 20 };

  typedef itk::Image<short, 3> ImageType;
  typedef itk::Image<short, 3> ProbabilityMapType;
  typedef itk::Image<unsigned char, 3> MaskType;
  typedef itk::Point<float, 3> PointType;

  struct ComparePair {
    typedef std::pair<PointType, float> ValueType;

    bool operator()(const ValueType &a, const ValueType &b) const {
      return a.second > b.second;
    }
  };

  unsigned int GetMaxDetections() const {
    return (1 << 31);
  }

  float GetDistanceThreshold() const {
    return 10.0f;
  }

  void SetDataFolder(const std::string &strDataFolder) {
    m_clProstateCAD.SetDataFolder(strDataFolder);
  }

  void SetProbabilityMapFolder(const std::string &strProbabilityMapFolder) {
    m_strProbabilityMapFolder = strProbabilityMapFolder;
  }

  bool LoadList(const std::string &strListFile);
  bool LoadNormalList(const std::string &strNormalListFile);

  std::string GetDataFolder() const {
    return m_clProstateCAD.GetDataFolder();
  }

  const std::string & GetProbabilityMapFolder() const {
    return m_strProbabilityMapFolder;
  }

  const std::vector<std::string> & GetPatientIDs() const {
    return m_vPatientIDs;
  }

  virtual ImageType::Pointer LoadT2WIVolume(const std::string &strPatientID) const {
    const std::string strVolumeFile = m_clProstateCAD.GetT2WIDataFolder() + '/' + strPatientID + ".%d.dcm";
    return nih::LoadDicomVolumeByIndex(strVolumeFile.c_str());
  }

  virtual ProbabilityMapType::Pointer LoadProbabilityMap(const std::string &strPatientID) const {
    const std::string strProbMapFile = GetProbabilityMapFolder() + '/' + strPatientID + ".mhd";
    return nih::LoadVolume<short>(strProbMapFile.c_str());
  }

  virtual MaskType::Pointer LoadProstate(const std::string &strPatientID);
  virtual MaskType::Pointer LoadCentralGland(const std::string &strPatientID);

  virtual bool Run();
  virtual bool Run(const std::string &strPatientID);
  virtual bool RunNormal(const std::string &strPatientID);

  bool NonMaximumSuppression(const std::string &strPatientID, std::vector<std::pair<PointType, float> > &vDetectionPairs);

  bool LoadBiopsyData(const std::string &strPatientID, std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints, const ImageType::SpacingType &clSpacing);
  int ComputeNearestNeighbor(float x, float y, float z, const std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints) const;

  static float Distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return std::sqrt(std::pow(x0-x1, 2) + std::pow(y0-y1, 2) + std::pow(z0-z1, 2));
  }

  void PrintOverall();
  void PrintGleason(unsigned int uiGleason);
  void PrintGreaterThanGleason(unsigned int uiGleason);

private:
  nih::ProstateCAD m_clProstateCAD;
  std::string m_strProbabilityMapFolder;
  std::vector<std::string> m_vPatientIDs, m_vNormalPatientIDs;

  DetectionHistogram m_a_clTruePositiveHistogramTZ[11];
  DetectionHistogram m_a_clTruePositiveHistogramPZ[11];

  DetectionHistogram m_clFalsePositiveHistogramTZ;
  DetectionHistogram m_clFalsePositiveHistogramPZ;
};

int main(int argc, char **argv) {
  const char * const p_cArg0 = argv[0];

  std::string strDataFolder;
  std::string strBiopsyList;
  std::string strNormalList;
  std::string strProbMapFolder;

  int c = 0;
  while ((c = getopt(argc, argv, "d:hl:n:p:")) != -1) {
    switch (c) {
    case 'd':
      strDataFolder = optarg;
      break;
    case 'h':
      Usage(p_cArg0);
      break;
    case 'l':
      strBiopsyList = optarg;
      break;
    case 'n':
      strNormalList = optarg;
      break;
    case 'p':
      strProbMapFolder = optarg;
      break;
    case '?':
    default:
      Usage(p_cArg0);
      break;
    }
  }

  argc -= optind;
  argv += optind;

  if (strBiopsyList.empty() || strProbMapFolder.empty() || strDataFolder.empty() || strNormalList.empty()) {
    Usage(p_cArg0); // Exits
  }

  EvalGleasonFROC clEvaluate;

  clEvaluate.SetDataFolder(strDataFolder);
  clEvaluate.SetProbabilityMapFolder(strProbMapFolder);

  if (!clEvaluate.LoadList(strBiopsyList)) {
    std::cerr << "Error: Could not load list file '" << strBiopsyList << "'." << std::endl;
    return -1;
  }

  if (!clEvaluate.LoadNormalList(strNormalList)) {
    std::cerr << "Error: Could not load normal list file '" << strNormalList << "'." << std::endl;
    return -1;
  }

  clEvaluate.Run();

  return 0;
}

bool EvalGleasonFROC::LoadList(const std::string &strListFile) {
  m_vPatientIDs.clear();

  std::ifstream listStream(strListFile.c_str());

  if (!listStream) {
    std::cerr << "Error: Could not load '" << strListFile << "'." << std::endl;
    return false;
  }

  std::string strLine;
  while (std::getline(listStream, strLine)) {
    size_t p = strLine.find('\r');

    if (p != std::string::npos)
      strLine.erase(p);

    if (strLine.empty())
      continue;

    m_vPatientIDs.push_back(strLine);
  }

  return m_vPatientIDs.size() > 0;
}

bool EvalGleasonFROC::LoadNormalList(const std::string &strListFile) {
  m_vNormalPatientIDs.clear();

  std::ifstream listStream(strListFile.c_str());

  if (!listStream) {
    std::cerr << "Error: Could not load '" << strListFile << "'." << std::endl;
    return false;
  }

  std::string strLine;
  while (std::getline(listStream, strLine)) {
    size_t p = strLine.find('\r');

    if (p != std::string::npos)
      strLine.erase(p);

    if (strLine.empty())
      continue;

    m_vNormalPatientIDs.push_back(strLine);
  }

  return m_vNormalPatientIDs.size() > 0;
}

EvalGleasonFROC::MaskType::Pointer EvalGleasonFROC::LoadProstate(const std::string &strPatientID) {
  ImageType::Pointer p_clImage = LoadT2WIVolume(strPatientID);

  if (!p_clImage)
    return MaskType::Pointer();

  MaskType::Pointer p_clMask = MaskType::New();

  p_clMask->SetOrigin(p_clImage->GetOrigin());
  p_clMask->SetSpacing(p_clImage->GetSpacing());
  p_clMask->SetRegions(p_clImage->GetBufferedRegion());

  p_clMask->Allocate(true);

  std::vector<unsigned int> vSlices; // Unused
  if (!m_clProstateCAD.LoadProstate(strPatientID, p_clMask, vSlices))
    return MaskType::Pointer();

  p_clMask->SetSpacing(p_clImage->GetSpacing());
  p_clMask->SetRegions(p_clImage->GetBufferedRegion());

  return p_clMask;
}

EvalGleasonFROC::MaskType::Pointer EvalGleasonFROC::LoadCentralGland(const std::string &strPatientID) {
  ImageType::Pointer p_clImage = LoadT2WIVolume(strPatientID);

  if (!p_clImage)
    return MaskType::Pointer();

  MaskType::Pointer p_clMask = MaskType::New();

  p_clMask->SetOrigin(p_clImage->GetOrigin());
  p_clMask->SetSpacing(p_clImage->GetSpacing());
  p_clMask->SetRegions(p_clImage->GetBufferedRegion());

  p_clMask->Allocate(true);

  std::vector<unsigned int> vSlices; // Unused
  if (!m_clProstateCAD.LoadCentralGland(strPatientID, p_clMask, vSlices))
    return MaskType::Pointer();

  p_clMask->SetSpacing(p_clImage->GetSpacing());
  p_clMask->SetRegions(p_clImage->GetBufferedRegion());

  return p_clMask;
}

bool EvalGleasonFROC::Run() {
  if (m_vPatientIDs.empty() || m_vNormalPatientIDs.empty()) {
    std::cerr << "Error: Biopsy and/or normal patients are missing.\n" << std::endl;
    return false;
  }

  for (int i = 0; i < 11; ++i) {
    m_a_clTruePositiveHistogramTZ[i].Clear();
    m_a_clTruePositiveHistogramPZ[i].Clear();
  }

  m_clFalsePositiveHistogramTZ.Clear();
  m_clFalsePositiveHistogramPZ.Clear();

  std::cout << "\nInfo: Processing biopsy cases ...\n" << std::endl;

  unsigned int uiPatientCount = 0;
  for (size_t i = 0; i < m_vPatientIDs.size(); ++i) {
    const std::string &strPatientID = m_vPatientIDs[i];

    std::cout << "Info: Processing '" << strPatientID << "' ..." << std::endl;

    if (!Run(strPatientID)) {
      std::cerr << "Error: Failed to process patient." << std::endl;
      continue;
    }

    ++uiPatientCount;
  }

  std::cout << "\nInfo: Processing normal cases ...\n" << std::endl;

#if 0
  uiPatientCount = 0;
  for (size_t i = 0; i < m_vNormalPatientIDs.size(); ++i) {
    const std::string &strPatientID = m_vNormalPatientIDs[i];

    std::cout << "Info: Processing '" << strPatientID << "' ..." << std::endl;

    if (!RunNormal(strPatientID)) {
      std::cerr << "Error: Failed to process patient." << std::endl;
      continue;
    }

    ++uiPatientCount;
  }
#endif

  m_clFalsePositiveHistogramTZ.SetPatientCount(uiPatientCount);
  m_clFalsePositiveHistogramPZ.SetPatientCount(uiPatientCount);

  std::cout << "Info: Performance:\n" << std::endl;
  PrintOverall();

  return true;
}

bool EvalGleasonFROC::Run(const std::string &strPatientID) {
  std::vector<nih::ProstateCAD::BiopsyPointType> vBiopsyPoints;

  ImageType::Pointer p_clT2WIVol = LoadT2WIVolume(strPatientID);
  if (!p_clT2WIVol) {
    std::cerr << "Error: Could not load T2WI volume for patient '" << strPatientID << "'." << std::endl;
    return false;
  }

  const ImageType::SpacingType &clSpacing = p_clT2WIVol->GetSpacing();

  if (!LoadBiopsyData(strPatientID, vBiopsyPoints, clSpacing)) {
    std::cerr << "Error: Could not load biopsy points for patient '" << strPatientID << "'." << std::endl;
    return false;
  }

  MaskType::Pointer p_clProstate = LoadProstate(strPatientID);
  if (!p_clProstate) {
    std::cerr << "Error: Failed to load prostate." << std::endl;
    return false;
  }

  MaskType::Pointer p_clCentralGland = LoadCentralGland(strPatientID);
  if (!p_clCentralGland) {
    std::cerr << "Error: Failed to load central gland." << std::endl;
    return false;
  }

  {
    std::vector<nih::ProstateCAD::BiopsyPointType> vTmpPoints;

    // Remove bad biopsy points (outside of prostate)
    for (size_t i = 0; i < vBiopsyPoints.size(); ++i) {
      const nih::ProstateCAD::BiopsyPointType &clBiopsyPoint = vBiopsyPoints[i];

      const int x = (int)(clBiopsyPoint.x / clSpacing[0] + 0.5f);
      const int y = (int)(clBiopsyPoint.y / clSpacing[1] + 0.5f);
      const int z = (int)(clBiopsyPoint.z / clSpacing[2] + 0.5f);

      const MaskType::IndexType clIndex = { x, y, z };

      if (p_clProstate->GetPixel(clIndex) != 0)
        vTmpPoints.push_back(clBiopsyPoint);
      else
        std::cout << "Info: Removing biopsy point: " << clBiopsyPoint.x << ", " << clBiopsyPoint.y << ", " << clBiopsyPoint.z << std::endl;
    }

    vBiopsyPoints.swap(vTmpPoints);
  }

  if (vBiopsyPoints.empty()) {
    std::cerr << "Error: No remaining biopsy points." << std::endl;
    return false;
  }

  std::vector<std::pair<PointType, float> > vDetections;
  if (!NonMaximumSuppression(strPatientID, vDetections)) {
    std::cerr << "Error: Non-maximum suppression failed." << std::endl;
    return false;
  }

  //if (vDetections.size() > GetMaxDetections())
    //vDetections.resize(GetMaxDetections()); // Keep top detections

  DetectionGraph clGraph(vBiopsyPoints, vDetections);

  std::vector<int> vPairs;
  const float fProbSum = clGraph.Optimize(vPairs);
  std::cout << "Info: Optimal probability sum = " << fProbSum << std::endl;

  for (size_t i = 0; i < vPairs.size(); ++i) {
    const int iDetIndex = vPairs[i];
    const nih::ProstateCAD::BiopsyPointType &clBiopsyPoint = vBiopsyPoints[i];

    const int x = (int)(clBiopsyPoint.x / clSpacing[0] + 0.5f);
    const int y = (int)(clBiopsyPoint.y / clSpacing[1] + 0.5f);
    const int z = (int)(clBiopsyPoint.z / clSpacing[2] + 0.5f);
    const MaskType::IndexType clIndex = { x, y, z };

    const bool bTransitionZone = (p_clCentralGland->GetPixel(clIndex) != 0);

    DetectionHistogram &clTruePositiveHistogram = bTransitionZone ? m_a_clTruePositiveHistogramTZ[clBiopsyPoint.gleason] : m_a_clTruePositiveHistogramPZ[clBiopsyPoint.gleason];
    DetectionHistogram &clFalsePositiveHistogram = bTransitionZone ? m_clFalsePositiveHistogramTZ : m_clFalsePositiveHistogramPZ;

    if (iDetIndex < 0) {
      // Not detected
      clTruePositiveHistogram.Update(0.0f);
      std::cout << "Info: True positive " << i << ": Zone=" << (bTransitionZone ? "transition" : "peripheral") << " Gleason=" << clBiopsyPoint.gleason << " Probability=" << 0.0f << std::endl;
    } 
    else {
      const std::pair<PointType, float> &clPair = vDetections[iDetIndex];
      const float fProb = clPair.second;

      clTruePositiveHistogram.Update(fProb);
      std::cout << "Info: True positive " << i << ": Zone=" << (bTransitionZone ? "transition" : "peripheral") << " Gleason=" << clBiopsyPoint.gleason << " Probability=" << fProb << std::endl;
    }
  }

  for (size_t i = 0; i < vDetections.size(); ++i) {
    if (std::find(vPairs.begin(), vPairs.end(), (int)i) != vPairs.end())
      continue;

    const std::pair<PointType, float> &clPair = vDetections[i];
    const PointType &clPoint = clPair.first;
    const float fProb = clPair.second;

    float fDist = GetDistanceThreshold() + 1.0f;
    {
      const int iNearest = ComputeNearestNeighbor(clPoint[0], clPoint[1], clPoint[2], vBiopsyPoints);
      if (iNearest >= 0) {
        const nih::ProstateCAD::BiopsyPointType &clBiopsyPoint = vBiopsyPoints[iNearest];
        fDist = Distance(clPoint[0], clPoint[1], clPoint[2], clBiopsyPoint.x, clBiopsyPoint.y, clBiopsyPoint.z);
      }
    }

    const int x = (int)(clPoint[0] / clSpacing[0] + 0.5f);
    const int y = (int)(clPoint[1] / clSpacing[1] + 0.5f);
    const int z = (int)(clPoint[2] / clSpacing[2] + 0.5f);
    const MaskType::IndexType clIndex = { x, y, z };

    const bool bTransitionZone = (p_clCentralGland->GetPixel(clIndex) != 0);
    DetectionHistogram &clFalsePositiveHistogram = bTransitionZone ? m_clFalsePositiveHistogramTZ : m_clFalsePositiveHistogramPZ;

    if (fDist > GetDistanceThreshold()) {
      clFalsePositiveHistogram.Update(fProb);
      std::cout << "Info: False positive: Zone=" << (bTransitionZone ? "transition" : "peripheral") << " Probability=" << fProb << std::endl;
    }
  }

  return true;
}

bool EvalGleasonFROC::RunNormal(const std::string &strPatientID) {
  ImageType::Pointer p_clT2WIVol = LoadT2WIVolume(strPatientID);
  if (!p_clT2WIVol) {
    std::cerr << "Error: Could not load T2WI volume for patient '" << strPatientID << "'." << std::endl;
    return false;
  }

  const ImageType::SpacingType &clSpacing = p_clT2WIVol->GetSpacing();

  MaskType::Pointer p_clProstate = LoadProstate(strPatientID);
  if (!p_clProstate) {
    std::cerr << "Error: Failed to load prostate." << std::endl;
    return false;
  }

  MaskType::Pointer p_clCentralGland = LoadProstate(strPatientID);
  if (!p_clCentralGland) {
    std::cerr << "Error: Failed to load central gland." << std::endl;
    return false;
  }

  std::vector<std::pair<PointType, float> > vDetections;
  if (!NonMaximumSuppression(strPatientID, vDetections)) {
    std::cerr << "Error: Non-maximum suppression failed." << std::endl;
    return false;
  }

  std::cout << "Info: Number of detections: " << vDetections.size() << std::endl;

  for (size_t i = 0; i < vDetections.size() && i < GetMaxDetections(); ++i) {
    const PointType &clPoint = vDetections[i].first;
    const float fProb = vDetections[i].second;


    const int x = (int)(clPoint[0] / clSpacing[0] + 0.5f);
    const int y = (int)(clPoint[1] / clSpacing[1] + 0.5f);
    const int z = (int)(clPoint[2] / clSpacing[2] + 0.5f);
    const MaskType::IndexType clIndex = { x, y, z };

    const bool bTransitionZone = (p_clCentralGland->GetPixel(clIndex) != 0);

    DetectionHistogram &clFalsePositiveHistogram = bTransitionZone ? m_clFalsePositiveHistogramTZ : m_clFalsePositiveHistogramPZ;

    clFalsePositiveHistogram.Update(fProb);
  }

  return true;
}

bool EvalGleasonFROC::NonMaximumSuppression(const std::string &strPatientID, std::vector<std::pair<PointType, float> > &vDetectionPairs) {
  vDetectionPairs.clear();

  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);
  if (!p_clProbMap) {
    std::cerr << "Error: Could not load probability map for patient '" << strPatientID << "'." << std::endl;
    return false;
  }

  ImageType::Pointer p_clT2WIVol = LoadT2WIVolume(strPatientID);

  const ProbabilityMapType::SpacingType &clSpacing = p_clT2WIVol->GetSpacing();

  std::cout << "Info: Probability map spacing: " << clSpacing << std::endl;

  const unsigned int uiWindowWidth = (unsigned int)(2*GetDistanceThreshold() / clSpacing[0] + 0.5f);
  const unsigned int uiWindowHeight = (unsigned int)(2*GetDistanceThreshold() / clSpacing[1] + 0.5f);
  const unsigned int uiWindowDepth = (unsigned int)(2*GetDistanceThreshold() / clSpacing[2] + 0.5f);

  ProbabilityMapType::Pointer p_clDetMap = nih::NonMaximumSuppression3D<short>(p_clProbMap, uiWindowWidth, uiWindowHeight, uiWindowDepth);
  if (!p_clDetMap) {
    std::cerr << "Error: NMS failed." << std::endl;
    return false;
  }

  const ProbabilityMapType::SizeType &clSize = p_clDetMap->GetBufferedRegion().GetSize();

  for (int z = 0; z < clSize[2]; ++z) {
    const float fZ = clSpacing[2]*z;

    for (int y = 0; y < clSize[1]; ++y) {
      const float fY = clSpacing[1]*y;

      for (int x = 0; x < clSize[0]; ++x) {
        const float fX = clSpacing[0]*x;
        const ProbabilityMapType::IndexType clIndex = { x, y, z };

        const short sProb = p_clDetMap->GetPixel(clIndex);

        if (sProb <= 0)
          continue;

        float fProb = sProb / 1024.0f;

        PointType clPoint;
        clPoint[0] = fX;
        clPoint[1] = fY;
        clPoint[2] = fZ;

#if 0
        // Compute momentum in some way

        const int xBegin = std::max(0, (int)(x - uiWindowWidth/2));
        const int xEnd = std::min((int)(clSize[0]), (int)(xBegin + uiWindowWidth));

        const int yBegin = std::max(0, (int)(y - uiWindowHeight/2));
        const int yEnd = std::min((int)(clSize[1]), (int)(yBegin + uiWindowHeight));

        const int zBegin = std::max(0, (int)(z - uiWindowDepth/2));
        const int zEnd = std::min((int)(clSize[2]), (int)(zBegin + uiWindowDepth));

        std::vector<float> vProbs;
        //fProb = 0.0f;
        unsigned int uiCount = 0;

        for (int k = zBegin; k < zEnd; ++k) {
          for (int j = yBegin; j < yEnd; ++j) {
            for (int i = xBegin; i < xEnd; ++i) {
              const ProbabilityMapType::IndexType clWindowIndex = { i, j, k };
              const short sWindowProb = p_clProbMap->GetPixel(clWindowIndex);

              if (sWindowProb <= 0)
                continue;

              const float fWindowProb = sWindowProb / 1024.0f;

              vProbs.push_back(fWindowProb);
              //fProb += fWindowProb;
              ++uiCount;
            }
          }
        }

        const float fPerc = 0.9f;

        std::vector<float>::iterator medItr = vProbs.begin() + (size_t)(fPerc*vProbs.size());
        std::nth_element(vProbs.begin(), medItr, vProbs.end());

        fProb = *medItr;

        if (std::distance(medItr, vProbs.end()) < (3*3*3)/(clSpacing[0]*clSpacing[1]*clSpacing[2]))
          continue;

        //if (uiCount > 0)
          //fProb /= uiCount;
#endif

        if (fProb > 0.0f)
          vDetectionPairs.push_back(std::make_pair(clPoint, fProb));

      }
    }
  }

  std::sort(vDetectionPairs.begin(), vDetectionPairs.end(), ComparePair());

#if 0
  for (size_t i = 0; i < vDetectionPairs.size(); ++i) {
    const PointType &clPoint = vDetectionPairs[i].first;
    float &fProb = vDetectionPairs[i].second;

    const int x = (int)(clPoint[0] / clSpacing[0]);
    const int y = (int)(clPoint[1] / clSpacing[1]);
    const int z = (int)(clPoint[2] / clSpacing[2]);
    const ProbabilityMapType::IndexType clIndex = { x, y, z };

    fProb = p_clProbMap->GetPixel(clIndex) / 1024.0f;
  }
#endif


  return true;
}

bool EvalGleasonFROC::LoadBiopsyData(const std::string &strPatientID, std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints, const ImageType::SpacingType &clSpacing) {
  nih::ProstateCAD clCAD;

  vBiopsyPoints.clear();

  clCAD.SetDataFolder(GetDataFolder());

  std::vector<nih::ProstateCAD::BiopsyPointType> vTmpPoints;
  if (!clCAD.LoadBiopsyData(vTmpPoints, clSpacing))
    return false;

  for (size_t i = 0; i < vTmpPoints.size(); ++i) {
    if (vTmpPoints[i].name == strPatientID)
      vBiopsyPoints.push_back(vTmpPoints[i]);
  }

  return vBiopsyPoints.size() > 0;
}

int EvalGleasonFROC::ComputeNearestNeighbor(float x, float y, float z, const std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints) const {

  if (vBiopsyPoints.empty())
    return -1;

  int iMinIndex = 0;
  float fMinDistance = Distance(x,y,z,vBiopsyPoints[0].x, vBiopsyPoints[0].y, vBiopsyPoints[0].z);

  for (size_t i = 1; i < vBiopsyPoints.size(); ++i) {
    const float fDistance = Distance(x,y,z,vBiopsyPoints[i].x, vBiopsyPoints[i].y, vBiopsyPoints[i].z);

    if (fDistance < fMinDistance) {
      iMinIndex = (int)i;
      fMinDistance = fDistance;
    }
  }

  return iMinIndex;
}

void EvalGleasonFROC::PrintOverall() {

  std::cout << "\nInfo: Gleason-grade analysis:\n" << std::endl;

  PrintGleason(0);
  for (unsigned int i = 6; i <= 10; ++i)
    PrintGleason(i);

  std::cout << "\nInfo: Cumulative Gleason-grade analysis:\n" << std::endl;

  PrintGreaterThanGleason(0);
  for (unsigned int i = 6; i <= 10; ++i)
    PrintGreaterThanGleason(i);
}

void EvalGleasonFROC::PrintGleason(unsigned int uiGleason) {
  if (uiGleason > 10)
    return;

  std::cout << "\nInfo: FROC curve for Gleason grade " << uiGleason << ":\n" << std::endl;

  DetectionHistogram clOverallTPHistogram = m_a_clTruePositiveHistogramPZ[uiGleason];
  clOverallTPHistogram.Merge(m_a_clTruePositiveHistogramTZ[uiGleason]);

  DetectionHistogram clOverallFPHistogram = m_clFalsePositiveHistogramPZ;
  clOverallFPHistogram.Merge(m_clFalsePositiveHistogramTZ);

  clOverallFPHistogram.SetPatientCount(m_clFalsePositiveHistogramPZ.GetPatientCount());

  std::cout << "\nInfo: Peripheral Zone:" << std::endl;
  std::cout << BindFROC(m_a_clTruePositiveHistogramPZ[uiGleason], m_clFalsePositiveHistogramPZ) << std::endl;

  std::cout << "\nInfo: Total number of true positives: " << m_a_clTruePositiveHistogramPZ[uiGleason].GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << m_clFalsePositiveHistogramPZ.GetTotalCount() << ", " << m_clFalsePositiveHistogramPZ.GetPatientCount() << std::endl;

  std::cout << "\nInfo: Transition Zone:" << std::endl;
  std::cout << BindFROC(m_a_clTruePositiveHistogramTZ[uiGleason], m_clFalsePositiveHistogramTZ) << std::endl;

  std::cout << "\nInfo: Total number of true positives: " << m_a_clTruePositiveHistogramTZ[uiGleason].GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << m_clFalsePositiveHistogramTZ.GetTotalCount() << ", " << m_clFalsePositiveHistogramTZ.GetPatientCount() << std::endl;

  std::cout << "\nInfo: Overall:" << std::endl;
  std::cout << BindFROC(clOverallTPHistogram, clOverallFPHistogram) << std::endl;

  std::cout << "\nInfo: Total number: " << clOverallTPHistogram.GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << clOverallFPHistogram.GetTotalCount() << ", " << clOverallFPHistogram.GetPatientCount() << std::endl;

}

void EvalGleasonFROC::PrintGreaterThanGleason(unsigned int uiGleason) {
  if (uiGleason > 10)
    return;

  std::cout << "\nInfo: FROC curve for Gleason grades >= " << uiGleason << ":\n" << std::endl;

  DetectionHistogram clCumTPHistogramTZ, clCumTPHistogramPZ, clOverallCumTPHistogram;

  for (unsigned int i = uiGleason; i <= 10; ++i) {
    clCumTPHistogramTZ.Merge(m_a_clTruePositiveHistogramTZ[i]);
    clCumTPHistogramPZ.Merge(m_a_clTruePositiveHistogramPZ[i]);
  }

  clOverallCumTPHistogram = clCumTPHistogramPZ;
  clOverallCumTPHistogram.Merge(clCumTPHistogramTZ);

  DetectionHistogram clOverallFPHistogram = m_clFalsePositiveHistogramPZ;
  clOverallFPHistogram.Merge(m_clFalsePositiveHistogramTZ);

  clOverallFPHistogram.SetPatientCount(m_clFalsePositiveHistogramPZ.GetPatientCount());

  std::cout << "\nInfo: Peripheral Zone:" << std::endl;
  std::cout << BindFROC(clCumTPHistogramPZ, m_clFalsePositiveHistogramPZ) << std::endl;

  std::cout << "\nInfo: Total number of true positives: " << clCumTPHistogramPZ.GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << m_clFalsePositiveHistogramPZ.GetTotalCount() << ", " << m_clFalsePositiveHistogramPZ.GetPatientCount() << std::endl;

  std::cout << "\nInfo: Transition Zone:" << std::endl;
  std::cout << BindFROC(clCumTPHistogramTZ, m_clFalsePositiveHistogramTZ) << std::endl;

  std::cout << "\nInfo: Total number of true positives: " << clCumTPHistogramTZ.GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << m_clFalsePositiveHistogramTZ.GetTotalCount() << ", " << m_clFalsePositiveHistogramTZ.GetPatientCount() << std::endl;

  std::cout << "\nInfo: Overall:" << std::endl;
  std::cout << BindFROC(clOverallCumTPHistogram, clOverallFPHistogram) << std::endl;

  std::cout << "\nInfo: Total number: " << clOverallCumTPHistogram.GetTotalCount() << std::endl;
  std::cout << "Info: Total number of false positives, patients: " << clOverallFPHistogram.GetTotalCount() << ", " << clOverallFPHistogram.GetPatientCount() << std::endl;
}
  

float DetectionGraph::OptimizeHelper(int iBiopsyPoint, float fProbSum, std::vector<int> &vPairs, std::vector<int> &vClaimed) const {
  if (iBiopsyPoint < 0 || iBiopsyPoint >= m_vBiopsyPoints.size())
    return fProbSum;

  float fMaxProbSum = fProbSum;

  const int iGleason = m_vBiopsyPoints[iBiopsyPoint].gleason;
  //const float fWeight = 1.0f + std::exp(0.41f*(iGleason - 6));
  const float fWeight = std::pow(1.5f, iGleason-6);

  std::vector<int> vBestPairs(vPairs);

  vPairs[iBiopsyPoint] = -1;

  // Try with no connection
  float fNewProbSum = OptimizeHelper(iBiopsyPoint+1, fProbSum, vPairs, vClaimed);

  if (fNewProbSum > fMaxProbSum) {

    fMaxProbSum = fNewProbSum;
    vBestPairs = vPairs;
  }

  const std::vector<int> &vNeighbors = m_vGraph[iBiopsyPoint];

  for (size_t i = 0; i < vNeighbors.size(); ++i) {
    const int iDetIndex = vNeighbors[i];

    if (vClaimed[iDetIndex])
      continue;

    const float fDetProb = fWeight * m_vDetectionPairs[iDetIndex].second;
    //const float fDetProb = m_vDetectionPairs[iDetIndex].second;
    //const float fDetProb = 1.0f;

    vClaimed[iDetIndex] = 1;
    vPairs[iBiopsyPoint] = iDetIndex;
    fNewProbSum = OptimizeHelper(iBiopsyPoint+1, fProbSum + fDetProb, vPairs, vClaimed);
    vClaimed[iDetIndex] = 0;



    if (fNewProbSum > fMaxProbSum) {
      fMaxProbSum = fNewProbSum;
      vBestPairs = vPairs;
    }
  }

  vPairs = vBestPairs;

  return fMaxProbSum;
}

void DetectionGraph::ComputeGraph() {
  m_vGraph.clear();

  m_vGraph.resize(m_vBiopsyPoints.size());

  for (size_t i = 0; i < m_vGraph.size(); ++i) {
    const nih::ProstateCAD::BiopsyPointType &clBiopsyPoint = m_vBiopsyPoints[i];
    std::vector<int> &vNeighbors = m_vGraph[i];

    for (size_t j = 0; j < m_vDetectionPairs.size(); ++j) {
      const PointType &clPoint = m_vDetectionPairs[j].first;
      const float fDistance = Distance(clPoint[0], clPoint[1], clPoint[2], clBiopsyPoint.x, clBiopsyPoint.y, clBiopsyPoint.z);

      if (fDistance <= DistanceThreshold())
        vNeighbors.push_back((int)j);
    }
  }
}
