#include <cstdlib>
#include <cmath>
#include <sstream>
#include <fstream>
#include <memory>
#include "ROCCurve.h"
#include "bsdgetopt.h"
#include "ProstateCAD.h"
#include "Morphology.h"
#include "MeanShift.h"

class MuteStream {
public:
  MuteStream(std::ostream &os)
  : m_os(os) {
    m_p_oldBuf = os.rdbuf(m_blackHole.rdbuf());
  }

  ~MuteStream() {
    m_os.rdbuf(m_p_oldBuf);
  }

private:
  std::stringstream m_blackHole;
  std::streambuf *m_p_oldBuf;
  std::ostream &m_os;
};

class Evaluation {
public:
  typedef itk::Image<short, 3> ImageType;
  typedef itk::Image<short, 3> ProbabilityMapType;
  typedef itk::Image<unsigned char, 3> MaskType;
  typedef nih::ProstateCAD::BiopsyPointType BiopsyPointType;

  Evaluation() {
    m_bLabRocOutput = false;
  }

  virtual ~Evaluation() { }

  void SetDataFolder(const std::string &strDataFolder) {
    m_clProstateCAD.SetDataFolder(strDataFolder);
  }

  void SetProbabilityMapFolder(const std::string &strProbabilityMapFolder) {
    m_strProbabilityMapFolder = strProbabilityMapFolder;
  }

  bool LoadList(const std::string &strListFile);

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
  virtual MaskType::Pointer LoadTumors(const std::string &strPatientID);
  virtual std::vector<BiopsyPointType> LoadBiopsies(const std::string &strPatientID, const ImageType::SpacingType &clSpacing);

  virtual bool RunLabRoc();
  virtual bool Run();

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) = 0;
  virtual void PrintROC(nih::ROCCurve &clROC, unsigned int /*uiNumCases*/) const {
    std::cout << nih::ROCFormat(clROC) << std::endl;
    std::cout << "AUC: " << clROC.ComputeAUC() << std::endl;
  }

  virtual bool GetDetections(const std::string &strPatientId, std::vector<float> &vProbsNegatives, std::vector<float> &vProbsPositives) {
    std::cerr << "Error: GetDetections() is not implemented." << std::endl;
    return false;
  }

  virtual void SetLabRocOutput(bool bOutput) {
    m_bLabRocOutput = bOutput;
  }

  virtual bool GetLabRocOutput() const {
    return m_bLabRocOutput;
  }

private:
  nih::ProstateCAD m_clProstateCAD;
  std::string m_strProbabilityMapFolder;

  std::vector<std::string> m_vPatientIDs;
  bool m_bLabRocOutput;
};

class EvaluationFactory {
public:
  class Node {
  public:
    virtual ~Node() { }
    virtual std::shared_ptr<Evaluation> Create() = 0;
  };

  typedef std::map<std::string, std::shared_ptr<Node> > MapType;

  template<typename EvaluationType>
  class NodeTemplate : public Node {
  public:
    virtual ~NodeTemplate() { }
    virtual std::shared_ptr<Evaluation> Create() {
      return std::make_shared<EvaluationType>();
    }
  };

  static EvaluationFactory & GetInstance() {
    static EvaluationFactory clSingleton;
    return clSingleton;
  }

  template<typename EvaluationType>
  void Register(const std::string &strName) {
    m_mNodes[strName].reset(new NodeTemplate<EvaluationType>());
  }

  std::shared_ptr<Evaluation> Create(const std::string &strName) {
    MapType::iterator itr = m_mNodes.find(strName);
    return itr != m_mNodes.end() ? itr->second->Create() : std::shared_ptr<Evaluation>();
  }

private:
  MapType m_mNodes;

  EvaluationFactory() { }
};

class VoxelPerformanceEvaluation : public Evaluation {
public:
  virtual ~VoxelPerformanceEvaluation() { }

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC);
  virtual void PrintROC(nih::ROCCurve &clROC, unsigned int /*uiNumCases*/) const {
    std::cout << nih::PrecisionRecallFormat(clROC) << std::endl;
  }
};

class LitjensPerformanceEvaluation : public Evaluation {
public:

  template<typename FirstType, typename SecondType>
  struct ComparePair {
    typedef std::pair<FirstType, SecondType> ValueType;

    bool operator()(const ValueType &a, const ValueType &b) const {
      return a.second < b.second;
    }
  };

  typedef itk::Point<float, 3> PointType;

  virtual ~LitjensPerformanceEvaluation() { }

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC);

  virtual void PrintROC(nih::ROCCurve &clROC, unsigned int uiNumCases) const {
    std::cout << nih::FROCCurveFormat(clROC, uiNumCases) << std::endl;
  }

  virtual std::vector<PointType> LoadCancerPoints(const std::string &strPatientID);

  template<typename IteratorType>
  std::pair<IteratorType, float> FindNearest(IteratorType begin, IteratorType end, const PointType &clPoint);

  template<typename IteratorType>
  std::vector<std::pair<IteratorType, float> > RankNearest(IteratorType begin, IteratorType end, const PointType &clPoint);

  template<typename VoxelType>
  void ComputeCCCentroids(std::vector<PointType> &vCentroids, typename itk::Image<VoxelType, 3>::Pointer p_clMask) const;

  void FindExtrema(std::vector<PointType> &vCentroids, ProbabilityMapType::Pointer p_clProbMap, int iWidth, int iHeight, int iDepth) const;
};

// Mean shift
class TumorPerformanceEvaluation : public Evaluation {
public:
  typedef itk::Point<float, 3> PointType;

  virtual ~TumorPerformanceEvaluation() { }

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC);

  virtual bool GetDetections(const std::string &strPatientId, std::vector<float> &vProbsNegatives, std::vector<float> &vProbsPositives);
  //virtual void PrintROC(nih::ROCCurve &clROC, unsigned int uiNumCases) const {
    //std::cout << nih::FROCCurveFormat(clROC, uiNumCases) << std::endl;
  //}
};

class BiopsyPerformanceEvaluation : public Evaluation {
public:
  virtual ~BiopsyPerformanceEvaluation() { }

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC);

  virtual bool LoadBiopsyData(const std::string &strPatientID, std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints, const ImageType::SpacingType &clSpacing);

  int ComputeNearestNeighbor(float x, float y, float z, const std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints) const;

  static float Distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return std::sqrt(std::pow(x0-x1, 2) + std::pow(y0-y1, 2) + std::pow(z0-z1, 2));
  }
};

class NMSPerformanceEvaluation : public LitjensPerformanceEvaluation {
public:
  virtual ~NMSPerformanceEvaluation() { }

  virtual bool ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC);

private:
  static float GetDistanceThreshold() { return 10.0f; }

  void DoROC(std::vector<std::pair<PointType, float> > &vDetections, std::vector<PointType> &vTumors, nih::ROCCurve &clROC);
};

template<typename SuperClassT>
class EvaluationMasksOnly : public SuperClassT {
public:
  typedef SuperClassT SuperClass;
  typedef typename SuperClass::ImageType ImageType;
  typedef typename SuperClass::ProbabilityMapType ProbabilityMapType;
  typedef typename SuperClass::MaskType MaskType;

  using SuperClass::GetProbabilityMapFolder;
  using SuperClass::GetDataFolder;

  virtual ~EvaluationMasksOnly() { }

  virtual typename ImageType::Pointer LoadT2WIVolume(const std::string &strPatientID) const {
    const std::string strProbMapFile = GetProbabilityMapFolder() + '/' + strPatientID + ".mhd";
    return nih::LoadVolume<short>(strProbMapFile.c_str());
  }

  virtual typename ProbabilityMapType::Pointer LoadProbabilityMap(const std::string &strPatientID) const {
    const std::string strProbMapFile = GetProbabilityMapFolder() + '/' + strPatientID + ".mhd";
    return nih::LoadVolume<short>(strProbMapFile.c_str());
  }

  virtual typename MaskType::Pointer LoadProstate(const std::string &strPatientID) {
    const std::string strVolumeFile = GetDataFolder() + "/Prostates/" + strPatientID + ".mhd";
    return nih::LoadVolume<unsigned char>(strVolumeFile.c_str());
  }

  virtual typename MaskType::Pointer LoadTumors(const std::string &strPatientID) {
    const std::string strVolumeFile = GetDataFolder() + "/Tumors/" + strPatientID + ".mhd";
    return nih::LoadVolume<unsigned char>(strVolumeFile.c_str());
  }
};

typedef EvaluationMasksOnly<LitjensPerformanceEvaluation> LitjensPerformanceEvaluationMasksOnly;
typedef EvaluationMasksOnly<TumorPerformanceEvaluation> TumorPerformanceEvaluationMasksOnly;
typedef EvaluationMasksOnly<VoxelPerformanceEvaluation> VoxelPerformanceEvaluationMasksOnly;

void Usage(const char *p_cArg0) {
  std::cerr << "Usage: " << p_cArg0 << " [-hL] -d dataFolder -l testList -p probabilityMapFolder [-t tumor]" << std::endl;
  exit(1);
}

int main(int argc, char **argv) {
  const char * const p_cArg0 = argv[0];

  EvaluationFactory &clFactory = EvaluationFactory::GetInstance();

  clFactory.Register<VoxelPerformanceEvaluation>("voxel");
  clFactory.Register<TumorPerformanceEvaluation>("tumor");
  clFactory.Register<LitjensPerformanceEvaluation>("litjens");
  clFactory.Register<BiopsyPerformanceEvaluation>("biopsy");
  clFactory.Register<NMSPerformanceEvaluation>("nms");
  clFactory.Register<VoxelPerformanceEvaluationMasksOnly>("voxelMask");
  clFactory.Register<TumorPerformanceEvaluationMasksOnly>("tumorMask");

  std::string strDataFolder;
  std::string strListFile;
  std::string strProbMapFolder;
  std::string strType = "tumor";

  bool bLabRocOutput = false;

  int c = 0;
  while ((c = getopt(argc, argv, "d:hl:Lp:t:")) != -1) {
    switch (c) {
    case 't':
      strType = optarg;
      break;
    case 'd':
      strDataFolder = optarg;
      break;
    case 'l':
      strListFile = optarg;
      break;
    case 'L':
      bLabRocOutput = true;
      break;
    case 'h':
      Usage(p_cArg0); // Exits
      break;
    case 'p':
      strProbMapFolder = optarg;
      break;
    case '?':
    default:
      Usage(p_cArg0); // Exits
      break;
    }
  }

  std::cerr << "Info: Evaluating with '" << strType << "' driver ..." << std::endl;

  std::shared_ptr<Evaluation> p_clEvaluator = clFactory.Create(strType);
  if (!p_clEvaluator) {
    std::cerr << "Error: Could not create evaluator of type '" << strType << "'." << std::endl;
    return -1;
  }

  p_clEvaluator->SetDataFolder(strDataFolder);
  p_clEvaluator->SetProbabilityMapFolder(strProbMapFolder);
  p_clEvaluator->SetLabRocOutput(bLabRocOutput);

  if (!p_clEvaluator->LoadList(strListFile))
    return -1;

  p_clEvaluator->Run();

  return 0;
}

bool Evaluation::LoadList(const std::string &strListFile) {
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

Evaluation::MaskType::Pointer Evaluation::LoadProstate(const std::string &strPatientID) {
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

Evaluation::MaskType::Pointer Evaluation::LoadTumors(const std::string &strPatientID) {
  ImageType::Pointer p_clImage = LoadT2WIVolume(strPatientID);

  if (!p_clImage)
    return MaskType::Pointer();

  MaskType::Pointer p_clMask = MaskType::New();

  p_clMask->SetRegions(p_clImage->GetBufferedRegion());
  p_clMask->SetSpacing(p_clImage->GetSpacing());
  p_clMask->SetRegions(p_clImage->GetBufferedRegion());

  p_clMask->Allocate(true);

  std::vector<unsigned int> vSlices; // Unused
  if (!m_clProstateCAD.LoadTumors(strPatientID, p_clMask, vSlices))
    return MaskType::Pointer();

  p_clMask->SetOrigin(p_clImage->GetOrigin());
  p_clMask->SetSpacing(p_clImage->GetSpacing());

  return p_clMask;
}

std::vector<Evaluation::BiopsyPointType> Evaluation::LoadBiopsies(const std::string &strPatientID, const ImageType::SpacingType &clSpacing) {
  std::vector<BiopsyPointType> vBiopsyPoints;

  if (!m_clProstateCAD.LoadBiopsyData(vBiopsyPoints, clSpacing))
    return std::vector<BiopsyPointType>();

  std::vector<BiopsyPointType> vPatientBiopsyPoints;
  for (size_t i = 0; i < vBiopsyPoints.size(); ++i) {
    if (vBiopsyPoints[i].name == strPatientID)
      vPatientBiopsyPoints.push_back(vBiopsyPoints[i]);
  }

  return vPatientBiopsyPoints;
}

bool Evaluation::RunLabRoc() {
  const float fPerc = 0.1f;
  std::vector<float> vAllPositiveDetections, vAllNegativeDetections;
  std::vector<float> vTmpPositiveDetections, vTmpNegativeDetections;

  {
    MuteStream clMute(std::cout);

    for (size_t i = 0; i < m_vPatientIDs.size(); ++i) {
      const std::string &strPatientID = m_vPatientIDs[i];

      vTmpPositiveDetections.clear();
      vTmpNegativeDetections.clear();

      if (!GetDetections(strPatientID, vTmpNegativeDetections, vTmpPositiveDetections)) {
        std::cerr << "Error: Failed to process patient." << std::endl;
        continue;
      }

      vAllPositiveDetections.insert(vAllPositiveDetections.end(), vTmpPositiveDetections.begin(), vTmpPositiveDetections.end());
      vAllNegativeDetections.insert(vAllNegativeDetections.end(), vTmpNegativeDetections.begin(), vTmpNegativeDetections.end());
    }
  }

  std::cout << "EvalProstateCAD\n";
  std::cout << "KIT\n";
  std::cout << "\"Prostate\"\n";
  std::cout << "CLL\n";

  const size_t negSize = (size_t)(fPerc*vAllNegativeDetections.size() + 0.5f);
  //const size_t posSize = (size_t)(fPerc*vAllPositiveDetections.size() + 0.5f);

  std::random_shuffle(vAllNegativeDetections.begin(), vAllNegativeDetections.end());
  //std::random_shuffle(vAllPositiveDetections.begin(), vAllPositiveDetections.end());

  vAllNegativeDetections.resize(negSize);
  //vAllPositiveDetections.resize(posSize);

  for (size_t i = 0; i < vAllNegativeDetections.size(); ++i)
    std::cout << vAllNegativeDetections[i] << '\n';

  std::cout << "*\n";

  for (size_t i = 0; i < vAllPositiveDetections.size(); ++i)
    std::cout << vAllPositiveDetections[i] << '\n';

  std::cout << '*' << std::endl;
  return true;
}

bool Evaluation::Run() {
  if (GetLabRocOutput())
    return RunLabRoc();

  nih::ROCCurve clOverallROC;

  for (size_t i = 0; i < m_vPatientIDs.size(); ++i) {
    const std::string &strPatientID = m_vPatientIDs[i];

    std::cout << "Info: Processing '" << strPatientID << "' ..." << std::endl;

    nih::ROCCurve clPatientROC;

    if (!ComputeROC(strPatientID, clPatientROC)) {
      std::cerr << "Error: Failed to process patient." << std::endl;
      continue;
    }

    std::cout << "Info: Patient performance:\n" << std::endl;

    PrintROC(clPatientROC, 1);

    clOverallROC.MergeCounts(clPatientROC);
  }

  std::cout << "Info: Overall performance:\n" << std::endl;
  PrintROC(clOverallROC, (unsigned int)m_vPatientIDs.size());

  return true;
}

bool VoxelPerformanceEvaluation::ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) {
  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);
  MaskType::Pointer p_clTumorMask = LoadTumors(strPatientID);
  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clTumorMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate, tumor, and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const MaskType::SizeType &clTumorMaskSize = p_clTumorMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();

  if (clProstateMaskSize != clTumorMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and tumor mask." << std::endl;
    return false;
  }

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  for (int z = 0; z < clProstateMaskSize[2]; ++z) {
    for (int y = 0; y < clProstateMaskSize[1]; ++y) {
      for (int x = 0; x < clProstateMaskSize[0]; ++x) {
        const MaskType::IndexType clIndex = { x, y, z };

        if (!p_clProstateMask->GetPixel(clIndex))
          continue;

        const int iGTLabel = p_clTumorMask->GetPixel(clIndex);
        const short sProb = p_clProbMap->GetPixel(clIndex);
        const float fProb = sProb/1024.0f;

        clROC.Update(fProb, iGTLabel);
      }
    }
  }

  return true;
}

bool LitjensPerformanceEvaluation::ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) {
  const float fDistThreshold = 10.0f;

  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);
  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SpacingType &clSpacing = p_clProstateMask->GetSpacing();

  std::cout << "Info: Spacing = " << clSpacing << std::endl;

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  const int iDilateX = (int)(0.5f*fDistThreshold / clSpacing[0] + 0.5f);
  const int iDilateY= (int)(0.5f*fDistThreshold / clSpacing[1] + 0.5f);
  const int iDilateZ = (int)(0.5f*fDistThreshold / clSpacing[2] + 0.5f);

  std::cout << "iDilateX = " << iDilateX << 
    ", iDilateY = " << iDilateY << 
    ", iDilateZ = " << iDilateZ << std::endl;

  std::vector<PointType> vTumors, vDetections;

  vTumors = LoadCancerPoints(strPatientID);
  FindExtrema(vDetections, p_clProbMap, iDilateX, iDilateY, iDilateZ);

  std::cout << "Info: Ground truth tumors:" << std::endl;
  for (size_t i = 0; i < vTumors.size(); ++i)
    std::cout << "Info: Tumor " << (i+1) << ": " << vTumors[i] << std::endl;

  std::cout << "Info: Detections:" << std::endl;
  for (size_t i = 0; i < vDetections.size(); ++i) {
    const PointType &clPoint = vDetections[i];

    std::cout << "Info: Detection " << (i+1) << ": " << vDetections[i] << std::endl;
  }

  std::vector<int> vDetected(vTumors.size(), 0);
  std::vector<int> vIsFalsePositive(vDetections.size(), 1); // Default to true

  // Phase 1, compute detections
  for (size_t i = 0; i < vTumors.size(); ++i) {
    const PointType &clTumorCentroid = vTumors[i];
    std::vector<std::pair<std::vector<PointType>::iterator, float> > vPairs = RankNearest(vDetections.begin(), vDetections.end(), clTumorCentroid);

    for (size_t j = 0; j < vPairs.size() && vPairs[j].second <= fDistThreshold; ++j) {
      const size_t detectionIndex = vPairs[j].first - vDetections.begin();

      const PointType &clPoint = *vPairs[j].first;

      const int x = std::min((int)(clPoint[0] / clSpacing[0] + 0.5f), (int)clProbMapSize[0]-1);
      const int y = std::min((int)(clPoint[1] / clSpacing[1] + 0.5f), (int)clProbMapSize[1]-1);
      const int z = std::min((int)(clPoint[2] / clSpacing[2] + 0.5f), (int)clProbMapSize[2]-1);

      const ProbabilityMapType::IndexType clIndex = { x, y, z };

      const short sProb = p_clProbMap->GetPixel(clIndex);
      const float fProb = std::min(sProb / 1024.0f, 1.0f);

      // If vIsFalsePositive[detectionIndex] is true, then this point isn't yet considered a detection for another lesion
      // Otherwise, it already detects another lesion

      if (vIsFalsePositive[detectionIndex]) {
        vDetected[i] = 1;
        vIsFalsePositive[detectionIndex] = 0;
        clROC.Update(fProb, 1);
        break;
      }
    }
  }

  // Anything else is a false negative
  for (size_t i = 0; i < vTumors.size(); ++i) {
    if (!vDetected[i]) {
      std::cout << "Info: Missed " << vTumors[i] << std::endl;

      clROC.Update(0.0f, 1);
    }
  }

  // Phase 2, prune detections to give any remaining false positives
  for (size_t i = 0; i < vTumors.size(); ++i) {
    const PointType &clTumorCentroid = vTumors[i];
    std::vector<std::pair<std::vector<PointType>::iterator, float> > vPairs = RankNearest(vDetections.begin(), vDetections.end(), clTumorCentroid);

    for (size_t j = 0; j < vPairs.size() && vPairs[j].second <= fDistThreshold; ++j) {
      const size_t detectionIndex = vPairs[j].first - vDetections.begin();
      vIsFalsePositive[detectionIndex] = 0; // No longer a false positive
    }
  }

  // Anything else is a false positive
  for (size_t i = 0; i < vDetections.size(); ++i) {
    const size_t detectionIndex = i;
    const PointType &clPoint = vDetections[i];

    if (!vIsFalsePositive[detectionIndex]) // Not a false positive
      continue;

    const int x = std::min((int)(clPoint[0] / clSpacing[0] + 0.5f), (int)clProbMapSize[0]-1);
    const int y = std::min((int)(clPoint[1] / clSpacing[1] + 0.5f), (int)clProbMapSize[1]-1);
    const int z = std::min((int)(clPoint[2] / clSpacing[2] + 0.5f), (int)clProbMapSize[2]-1);

    const ProbabilityMapType::IndexType clIndex = { x, y, z };

    const short sProb = p_clProbMap->GetPixel(clIndex);
    const float fProb = std::min(sProb / 1024.0f, 1.0f);

    clROC.Update(fProb, 0);
  }

  return true;
}

std::vector<LitjensPerformanceEvaluation::PointType> LitjensPerformanceEvaluation::LoadCancerPoints(const std::string &strPatientID) {
  std::vector<PointType> vPoints;
  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);

  if (!p_clProstateMask)
    return vPoints;

  const MaskType::SizeType &clSize = p_clProstateMask->GetBufferedRegion().GetSize();

  MaskType::Pointer p_clTumorMask = LoadTumors(strPatientID);

  if (p_clTumorMask.IsNotNull()) {
    for (int z = 0; z < clSize[2]; ++z) {
      for (int y = 0; y < clSize[1]; ++y) {
        for (int x = 0; x < clSize[0]; ++x) {
          const MaskType::IndexType clIndex = { x, y, z };

          if (!p_clProstateMask->GetPixel(clIndex))
            p_clTumorMask->SetPixel(clIndex, 0);
        }
      }
    }

    ComputeCCCentroids<unsigned char>(vPoints, p_clTumorMask);

    if (vPoints.size() > 0)
      return vPoints;
  }

  std::cout << "Info: Failed to load contours. Attempting to load biopsy points..." << std::endl;

  const ImageType::SpacingType &clSpacing = p_clProstateMask->GetSpacing();

  std::cout << "Info: Spacing: " << clSpacing << std::endl;

  std::vector<BiopsyPointType> vBiopsyPoints = LoadBiopsies(strPatientID, clSpacing);

  std::cout << "Info: Loaded " << vBiopsyPoints.size() << " biopsy points." << std::endl;

  for (size_t i = 0; i < vBiopsyPoints.size(); ++i) {
    PointType clPoint;
    clPoint[0] = vBiopsyPoints[i].x;
    clPoint[1] = vBiopsyPoints[i].y;
    clPoint[2] = vBiopsyPoints[i].z;

    const int x = (int)(clPoint[0] / clSpacing[0]);
    const int y = (int)(clPoint[1] / clSpacing[1]);
    const int z = (int)(clPoint[2] / clSpacing[2]);

    const MaskType::IndexType clIndex = { x, y, z };

    if (x >= 0 && x < clSize[0] && y >= 0 && y < clSize[1] && z >= 0 && z < clSize[2] && p_clProstateMask->GetPixel(clIndex) != 0)
      vPoints.push_back(clPoint); 
  }

  return vPoints;
}

template<typename IteratorType>
std::pair<IteratorType, float> LitjensPerformanceEvaluation::FindNearest(IteratorType begin, IteratorType end, const PointType &clPoint) {
  if (begin == end)
    return std::make_pair(end, 0.0f);

  std::pair<IteratorType, float> clNearest;

  clNearest.first = begin;
  clNearest.second = clPoint.EuclideanDistanceTo(*clNearest.first);

  ++begin;

  for ( ; begin != end; ++begin) {
    const float fDist = clPoint.EuclideanDistanceTo(*begin);

    if (fDist < clNearest.second) {
      clNearest.first = begin;
      clNearest.second = fDist;
    }
  }

  return clNearest;
}

template<typename IteratorType>
std::vector<std::pair<IteratorType, float> > LitjensPerformanceEvaluation::RankNearest(IteratorType begin, IteratorType end, const PointType &clPoint) {
  typedef ComparePair<IteratorType, float> CompareType;

  std::vector<std::pair<IteratorType, float> > vPairs;

  for ( ; begin != end; ++begin) {
    std::pair<IteratorType, float> clPair;

    clPair.first = begin;
    clPair.second = clPoint.EuclideanDistanceTo(*begin);

    vPairs.push_back(clPair);
  }

  std::sort(vPairs.begin(), vPairs.end(), CompareType());

  return vPairs;
}

template<typename VoxelType>
void LitjensPerformanceEvaluation::ComputeCCCentroids(std::vector<PointType> &vCentroids, typename itk::Image<VoxelType, 3>::Pointer p_clMask) const {
  vCentroids.clear();

  typedef itk::Image<VoxelType, 3> ImageType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::IndexType IndexType;

  const SizeType &clSize = p_clMask->GetBufferedRegion().GetSize();
  const SpacingType &clSpacing = p_clMask->GetSpacing();

  unsigned int uiNumComponents = 0;
  MaskType::Pointer p_clCCMask = nih::LabelConnectedComponents3D<VoxelType, unsigned char>(p_clMask, uiNumComponents);

  if (!p_clCCMask)
    return;

  vCentroids.resize(uiNumComponents, PointType(0.0f));

  std::vector<unsigned int> vCounts(vCentroids.size(), 0);

  for (int z = 0; z < clSize[2]; ++z) {
    for (int y = 0; y < clSize[1]; ++y) {
      for (int x = 0; x < clSize[0]; ++x) {
        const IndexType clIndex = { x, y, z };

        const unsigned char ucLabel = p_clCCMask->GetPixel(clIndex);

        if (ucLabel == 0)
          continue;

        if (ucLabel-1 >= vCentroids.size()) {
          std::cerr << "Error: Label out-of-bounds." << std::endl;
          continue; // Spam the console so we know
        }

        ++vCounts[ucLabel-1];
        vCentroids[ucLabel-1][0] += clSpacing[0]*x;
        vCentroids[ucLabel-1][1] += clSpacing[1]*y;
        vCentroids[ucLabel-1][2] += clSpacing[2]*z;
      }
    }
  }

  for (size_t i = 0; i < vCentroids.size(); ++i) {
    if (vCounts[i] != 0) {
      vCentroids[i][0] /= (float)vCounts[i];
      vCentroids[i][1] /= (float)vCounts[i];
      vCentroids[i][2] /= (float)vCounts[i];
    }
  }
}

void LitjensPerformanceEvaluation::FindExtrema(std::vector<PointType> &vExtrema, ProbabilityMapType::Pointer p_clProbMap, int iWidth, int iHeight, int iDepth) const {
  vExtrema.clear();

  typedef itk::FlatStructuringElement<3> StructuringElementType;
  typedef StructuringElementType::RadiusType RadiusType;

  RadiusType clRadius = { iWidth, iHeight, iDepth };

  StructuringElementType clBoxWithHole;
  clBoxWithHole.SetDecomposable(false);
  clBoxWithHole.SetRadius(clRadius);

  for (size_t i = 0; i < clBoxWithHole.Size(); ++i)
    clBoxWithHole[i] = true;

  clBoxWithHole[clBoxWithHole.GetCenterNeighborhoodIndex()] = false;

  ProbabilityMapType::Pointer p_clExtremaMap = nih::Dilate<short, StructuringElementType, 3>(p_clProbMap, clBoxWithHole);

  if (!p_clExtremaMap) {
    std::cerr << "Error: Could not compute 3D morphological dilation." << std::endl;
    return;
  }

  const ImageType::SpacingType &clSpacing = p_clProbMap->GetSpacing();

  std::cout << "*** Spacing = " << clSpacing << std::endl;

  const ProbabilityMapType::SizeType &clSize = p_clExtremaMap->GetBufferedRegion().GetSize();

  for (int z = 0; z < clSize[2]; ++z) {
    for (int y = 0; y < clSize[1]; ++y) {
      for (int x = 0; x < clSize[0]; ++x) {
        const ProbabilityMapType::IndexType clIndex = { x, y, z };

        const short sProb = p_clProbMap->GetPixel(clIndex);
        const short sProbDilated = p_clExtremaMap->GetPixel(clIndex);

        PointType clPoint;
        clPoint[0] = clSpacing[0]*x;
        clPoint[1] = clSpacing[1]*y;
        clPoint[2] = clSpacing[2]*z;

        if (sProb > sProbDilated)
          vExtrema.push_back(clPoint);

        p_clExtremaMap->SetPixel(clIndex, (sProb > sProbDilated) ? sProb : (short)0);
      }
    }
  }

  nih::SaveVolume<short>(p_clExtremaMap, "Extrema.mhd");

  exit(0);

  //ComputeCCCentroids<short>(vExtrema, p_clExtremaMap);
}

bool TumorPerformanceEvaluation::ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) {
  const float fTumorPerc = 0.9f;
  const float fCellPerc = 0.9f;
  const float fCellWidthMM = 3.0f;
  const float fCellHeightMM = 3.0f;
  const short sProbThreshold = (short)(0.1f*1024);

  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);
  MaskType::Pointer p_clTumorMask = LoadTumors(strPatientID);
  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clTumorMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate, tumor, and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const MaskType::SizeType &clTumorMaskSize = p_clTumorMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();
  const MaskType::SpacingType &clTumorMaskSpacing = p_clTumorMask->GetSpacing();

  if (clProstateMaskSize != clTumorMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and tumor mask." << std::endl;
    return false;
  }

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  unsigned int uiNumTumorCCs = 0;
  MaskType::Pointer p_clTumorCCMask = nih::LabelConnectedComponents3D<unsigned char, unsigned char>(p_clTumorMask, uiNumTumorCCs);

  std::vector<std::vector<short> > vAllProbs(uiNumTumorCCs+1);

  // Collect the 90th percentile of lesions in every lesion
  
  for (int z = 0; z < clTumorMaskSize[2]; ++z) {
    for (int y = 0; y < clTumorMaskSize[1]; ++y) {
      for (int x = 0; x < clTumorMaskSize[0]; ++x) {
        const MaskType::IndexType clIndex = { x, y, z };
        const unsigned char ucLabel = p_clTumorCCMask->GetPixel(clIndex);

        if (ucLabel == 0)
          continue;

        const short sProb = p_clProbMap->GetPixel(clIndex);

        vAllProbs[ucLabel].push_back(sProb);
      }
    }
  }

  for (size_t i = 1; i < vAllProbs.size(); ++i) {
    std::sort(vAllProbs[i].begin(), vAllProbs[i].end());

    const size_t index = (size_t)(fTumorPerc * vAllProbs[i].size());
    const short sProb = vAllProbs[i][index];
    const float fProb = sProb / 1024.0f;

    clROC.Update(fProb, 1);
  }

  const int iCellWidth = std::max(1, (int)(fCellWidthMM / clTumorMaskSpacing[0]));
  const int iCellHeight = std::max(1, (int)(fCellHeightMM / clTumorMaskSpacing[1]));

  const int iNumCellsX = clTumorMaskSize[0] / iCellWidth;
  const int iNumCellsY = clTumorMaskSize[1] / iCellHeight;

  std::vector<short> vCellProbs;

  for (int z = 0; z < clTumorMaskSize[2]; ++z) {
    for (int yCell = 0; yCell < iNumCellsY; ++yCell) {
      const int yBegin = yCell * iCellHeight;
      const int yEnd = std::min((int)clTumorMaskSize[1], yBegin + iCellHeight);

      for (int xCell = 0; xCell < iNumCellsX; ++xCell) {
        const int xBegin = xCell * iCellWidth;
        const int xEnd = std::min((int)clTumorMaskSize[0], xBegin + iCellWidth);

        vCellProbs.clear();
        bool bInTumor = false;

        for (int y = yBegin; y < yEnd; ++y) {
          for (int x = xBegin; x < xEnd; ++x) {
            const MaskType::IndexType clIndex = { x, y, z };
            const bool bInProstate = (p_clProstateMask->GetPixel(clIndex) != 0);            

            if (!bInProstate)
              continue;

            if (p_clTumorMask->GetPixel(clIndex) != 0) {
              bInTumor = true;
              break;
            }

            const short sProb = p_clProbMap->GetPixel(clIndex);
            vCellProbs.push_back(sProb);
          }
        }

        if (bInTumor || vCellProbs.empty()) // Cell contained part of a tumor or cell occurred outside of the prostate
          continue;

        std::sort(vCellProbs.begin(), vCellProbs.end());
        const size_t index = (size_t)(fCellPerc * vCellProbs.size());

        const float fProb = vCellProbs[index] / 1024.0f;

        clROC.Update(fProb, 0);
      }
    }
  }

  return true;
}

bool TumorPerformanceEvaluation::GetDetections(const std::string &strPatientID, std::vector<float> &vProbsNegatives, std::vector<float> &vProbsPositives) {
  const float fTumorPerc = 0.9f;
  const float fCellPerc = 0.9f;
  const float fCellWidthMM = 3.0f;
  const float fCellHeightMM = 3.0f;
  const short sProbThreshold = (short)(0.1f*1024);

  vProbsNegatives.clear();
  vProbsPositives.clear();

  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);
  MaskType::Pointer p_clTumorMask = LoadTumors(strPatientID);
  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clTumorMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate, tumor, and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const MaskType::SizeType &clTumorMaskSize = p_clTumorMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();
  const MaskType::SpacingType &clTumorMaskSpacing = p_clTumorMask->GetSpacing();

  if (clProstateMaskSize != clTumorMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and tumor mask." << std::endl;
    return false;
  }

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  unsigned int uiNumTumorCCs = 0;
  MaskType::Pointer p_clTumorCCMask = nih::LabelConnectedComponents3D<unsigned char, unsigned char>(p_clTumorMask, uiNumTumorCCs);

  std::vector<std::vector<short> > vAllProbs(uiNumTumorCCs+1);

  // Collect the 90th percentile of lesions in every lesion
  
  for (int z = 0; z < clTumorMaskSize[2]; ++z) {
    for (int y = 0; y < clTumorMaskSize[1]; ++y) {
      for (int x = 0; x < clTumorMaskSize[0]; ++x) {
        const MaskType::IndexType clIndex = { x, y, z };
        const unsigned char ucLabel = p_clTumorCCMask->GetPixel(clIndex);

        if (ucLabel == 0)
          continue;

        const short sProb = p_clProbMap->GetPixel(clIndex);

        vAllProbs[ucLabel].push_back(sProb);
      }
    }
  }

  for (size_t i = 1; i < vAllProbs.size(); ++i) {
    std::sort(vAllProbs[i].begin(), vAllProbs[i].end());

    const size_t index = (size_t)(fTumorPerc * vAllProbs[i].size());
    const short sProb = vAllProbs[i][index];
    const float fProb = sProb / 1024.0f;

    vProbsPositives.push_back(fProb);
  }

  const int iCellWidth = std::max(1, (int)(fCellWidthMM / clTumorMaskSpacing[0]));
  const int iCellHeight = std::max(1, (int)(fCellHeightMM / clTumorMaskSpacing[1]));

  const int iNumCellsX = clTumorMaskSize[0] / iCellWidth;
  const int iNumCellsY = clTumorMaskSize[1] / iCellHeight;

  std::vector<short> vCellProbs;

  for (int z = 0; z < clTumorMaskSize[2]; ++z) {
    for (int yCell = 0; yCell < iNumCellsY; ++yCell) {
      const int yBegin = yCell * iCellHeight;
      const int yEnd = std::min((int)clTumorMaskSize[1], yBegin + iCellHeight);

      for (int xCell = 0; xCell < iNumCellsX; ++xCell) {
        const int xBegin = xCell * iCellWidth;
        const int xEnd = std::min((int)clTumorMaskSize[0], xBegin + iCellWidth);

        vCellProbs.clear();
        bool bInTumor = false;

        for (int y = yBegin; y < yEnd; ++y) {
          for (int x = xBegin; x < xEnd; ++x) {
            const MaskType::IndexType clIndex = { x, y, z };
            const bool bInProstate = (p_clProstateMask->GetPixel(clIndex) != 0);            

            if (!bInProstate)
              continue;

            if (p_clTumorMask->GetPixel(clIndex) != 0) {
              bInTumor = true;
              break;
            }

            const short sProb = p_clProbMap->GetPixel(clIndex);
            vCellProbs.push_back(sProb);
          }
        }

        if (bInTumor || vCellProbs.empty()) // Cell contained part of a tumor or cell occurred outside of the prostate
          continue;

        std::sort(vCellProbs.begin(), vCellProbs.end());
        const size_t index = (size_t)(fCellPerc * vCellProbs.size());

        const float fProb = vCellProbs[index] / 1024.0f;

        vProbsNegatives.push_back(fProb);
      }
    }
  }

  return true;
}

bool BiopsyPerformanceEvaluation::ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) {
  const float fMinDistance = 5.0f;
  const float fPerc = 0.9f;

  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);

  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate, tumor, and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  const MaskType::SpacingType &clSpacing = p_clProstateMask->GetSpacing();

  std::vector<nih::ProstateCAD::BiopsyPointType> vBiopsyPoints;

  if (!LoadBiopsyData(strPatientID, vBiopsyPoints, clSpacing))
    return false;

  std::vector<std::vector<float> > vBiopsyProbs(vBiopsyPoints.size());

  for (int z = 0; z < clProbMapSize[2]; ++z) {
    for (int y = 0; y < clProbMapSize[1]; ++y) {
      for (int x = 0; x < clProbMapSize[0]; ++x) {
        const MaskType::IndexType clIndex = { x, y, z };

        if (!p_clProstateMask->GetPixel(clIndex))
          continue;

        const float fXmm = x*clSpacing[0];
        const float fYmm = y*clSpacing[1];
        const float fZmm = z*clSpacing[2];

        const int iIndex = ComputeNearestNeighbor(fXmm, fYmm, fZmm, vBiopsyPoints);

        if (iIndex < 0)
          continue; // Uhh?

        const nih::ProstateCAD::BiopsyPointType &clBiopsy = vBiopsyPoints[iIndex];
        std::vector<float> &vProbs = vBiopsyProbs[iIndex];
 
        const float fDistance = Distance(fXmm, fYmm, fZmm, clBiopsy.x, clBiopsy.y, clBiopsy.z);
        const float fProb = p_clProbMap->GetPixel(clIndex) / 1024.0f;

        if (fDistance <= fMinDistance)
          vProbs.push_back(fProb);
      }
    }
  }

  for (size_t i = 0; i < vBiopsyPoints.size(); ++i) {
    const int iLabel = vBiopsyPoints[i].label;
    std::vector<float> &vProbs = vBiopsyProbs[i];

    float fProb = 0.0f;

    if (vProbs.size() > 0) {
      std::sort(vProbs.begin(), vProbs.end());
      fProb = vProbs[(int)(vProbs.size() * fPerc)];
    }

    clROC.Update(fProb, iLabel);
  }

  return true;
}

bool BiopsyPerformanceEvaluation::LoadBiopsyData(const std::string &strPatientID, std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints, const ImageType::SpacingType &clSpacing) {
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

int BiopsyPerformanceEvaluation::ComputeNearestNeighbor(float x, float y, float z, const std::vector<nih::ProstateCAD::BiopsyPointType> &vBiopsyPoints) const {

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

bool NMSPerformanceEvaluation::ComputeROC(const std::string &strPatientID, nih::ROCCurve &clROC) {
  const float fDistThreshold = GetDistanceThreshold();

  MaskType::Pointer p_clProstateMask = LoadProstate(strPatientID);
  ProbabilityMapType::Pointer p_clProbMap = LoadProbabilityMap(strPatientID);

  if (!p_clProstateMask || !p_clProbMap) {
    std::cerr << "Error: Could not load prostate and/or probability map." << std::endl;
    return false;
  }

  const MaskType::SpacingType &clSpacing = p_clProstateMask->GetSpacing();

  std::cout << "Info: Spacing = " << clSpacing << std::endl;

  const MaskType::SizeType &clProstateMaskSize = p_clProstateMask->GetBufferedRegion().GetSize();
  const ProbabilityMapType::SizeType &clProbMapSize = p_clProbMap->GetBufferedRegion().GetSize();

  if (clProbMapSize != clProstateMaskSize) {
    std::cerr << "Error: Dimension mismatch between prostate and probability map." << std::endl;
    return false;
  }

  const unsigned int uiWindowWidthX = (int)(2*fDistThreshold / clSpacing[0] + 0.5f);
  const unsigned int uiWindowWidthY = (int)(2*fDistThreshold / clSpacing[1] + 0.5f);
  const unsigned int uiWindowWidthZ = (int)(2*fDistThreshold / clSpacing[2] + 0.5f);

  std::cout << "Info: NMS Window (in pixels): " << uiWindowWidthX << 'x' << uiWindowWidthY << 'x' << uiWindowWidthZ << std::endl;

  ProbabilityMapType::Pointer p_clNMSMap = nih::NonMaximumSuppression3D<short>(p_clProbMap, uiWindowWidthX, uiWindowWidthY, uiWindowWidthZ);

  //nih::SaveVolume<short>(p_clNMSMap, "NMSMap.mhd");

  std::vector<PointType> vTumors;
  std::vector<std::pair<PointType, float> > vDetections;

  vTumors = LoadCancerPoints(strPatientID);

  if (!p_clNMSMap) {
    std::cerr << "Error: Failed to compute NMS." << std::endl;
    return false;
  }

  if (vTumors.empty())
    std::cout << "Info: No ground truth tumors. All detections will be considered false positives." << std::endl;

  for (int z = 0; z < clProbMapSize[2]; ++z) {
    for (int y = 0; y < clProbMapSize[1]; ++y) {
      for (int x = 0; x < clProbMapSize[0]; ++x) {
        const ProbabilityMapType::IndexType clIndex = { x, y, z };

        if (p_clProstateMask->GetPixel(clIndex) == 0)
          continue;

        const float fProb = p_clNMSMap->GetPixel(clIndex) / 1024.0f;

        if (fProb > 0.0f) {
          PointType clPoint;

          clPoint[0] = x*clSpacing[0];
          clPoint[1] = y*clSpacing[1];
          clPoint[2] = z*clSpacing[2];

          vDetections.push_back(std::make_pair(clPoint, fProb));
        }
      }
    }
  }

  std::cout << "Info: There are " << vDetections.size() << " detections." << std::endl;

  const unsigned int uiNumTries = 100;
  nih::ROCCurve clBestROC;
  float fMaxAUC = -1.0f;

  for (unsigned int uiTry = 0; uiTry < uiNumTries; ++uiTry) {
    nih::ROCCurve clTmpROC;

    std::random_shuffle(vDetections.begin(), vDetections.end());

    DoROC(vDetections, vTumors, clTmpROC);

    const float fTmpAUC = clTmpROC.ComputeAUC();

    //if (uiTmpNumDetected > uiNumDetected) {
    if (fTmpAUC > fMaxAUC) {
      fMaxAUC = fTmpAUC;
      clBestROC = clTmpROC;
    }
  }

  clROC.MergeCounts(clBestROC);

  return true;
}

void NMSPerformanceEvaluation::DoROC(std::vector<std::pair<PointType, float> > &vDetections, std::vector<PointType> &vTumors, nih::ROCCurve &clROC) {
  const float fDistThreshold = GetDistanceThreshold();

  std::vector<bool> vAlreadyDetected(vTumors.size(), false), vFalsePositives(vDetections.size(), true);

  // First mark lesions detected and mark corresponding detections as not false positive
  for (size_t i = 0; i < vDetections.size(); ++i) {
    const PointType &clPoint = vDetections[i].first;
    const float fProb = vDetections[i].second;

    if (!vFalsePositives[i])
      continue;

    std::vector<std::pair<std::vector<PointType>::iterator, float> > vRanked = RankNearest(vTumors.begin(), vTumors.end(), clPoint);

    for (size_t j = 0; j < vRanked.size() && vRanked[j].second < fDistThreshold; ++j) {
      const size_t lesionIndex = vRanked[j].first - vTumors.begin();

      if (!vAlreadyDetected[lesionIndex]) {
        clROC.Update(fProb, 1);
        vAlreadyDetected[lesionIndex] = true;
        vFalsePositives[i] = false;
        break;
      }
    }
  }

  // Mark remaining nearby detections as not false positives
  for (size_t i = 0; i < vDetections.size(); ++i) {
    const PointType &clPoint = vDetections[i].first;

    if (!vFalsePositives[i])
      continue;

    std::vector<std::pair<std::vector<PointType>::iterator, float> > vRanked = RankNearest(vTumors.begin(), vTumors.end(), clPoint);

    if (vRanked.size() > 0 && vRanked.front().second < fDistThreshold)
      vFalsePositives[i] = false;
  }

  // Count false negatives
  for (size_t i = 0; i < vTumors.size(); ++i) {
    if (!vAlreadyDetected[i])
      clROC.Update(0.0f, 1);
  }

  // Count false positives
  for (size_t i = 0; i < vDetections.size(); ++i) {
    const float fProb = vDetections[i].second;

    if (vFalsePositives[i])
      clROC.Update(fProb, 0);
  }
}
