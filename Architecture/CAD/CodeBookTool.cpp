#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include "VoiFileIO.h"
#include "VoiToMask.h"
#include "CodeBookTree.h"
#include "ParameterContainer.h"
#include "ParameterMap.h"
#include "Common.h"
#include "getopt.h"

#undef LoadImage

void Usage(const char *p_cArg0) {
  std::cerr << "Usage: " << p_cArg0 << " train|convert|invert|plot [-h] -s minSampleSize -d dimension -i inputPixelType -o outputPixelType -l listFile -m foregroundExt -e imageExt -f codeBook.clf -r imageRoot -k maskRoot inputImage outputImage" << std::endl;
  exit(1);
}

template<typename PixelType>
PixelType StepValue(const PixelType & /*minValue*/, const PixelType & /*maxValue*/) {
  return PixelType(1);
}

float StepValue(const float &fMinValue, const float &fMaxValue) {
  return 0.01f*(fMaxValue - fMinValue);
}

double StepValue(const double &dMinValue, const double &dMaxValue) {
  return 0.01*(dMaxValue - dMinValue);
}

template<typename PixelType>
const char * GetPixelType() { return "Unknown"; }

#define DECLARE_PIXEL_TYPE(PixelTypeT) \
template<> \
const char * GetPixelType<PixelTypeT>() { return #PixelTypeT ; }

typedef unsigned short ushort;
typedef unsigned char uchar;
typedef unsigned int uint;

DECLARE_PIXEL_TYPE(char)
DECLARE_PIXEL_TYPE(uchar)
DECLARE_PIXEL_TYPE(short)
DECLARE_PIXEL_TYPE(ushort)
DECLARE_PIXEL_TYPE(int)
DECLARE_PIXEL_TYPE(uint)
DECLARE_PIXEL_TYPE(float)
DECLARE_PIXEL_TYPE(double)

template<typename PixelType, unsigned int Dimension>
struct VoiLoader {
  typedef itk::Image<PixelType, Dimension> ImageType;

  static typename ImageType::Pointer Load(const char *p_cFileName, const itk::Size<Dimension> &clSize) {
    std::cerr << "Error: Not implemented for pixel type " << GetPixelType<PixelType>() << " with dimension " << Dimension << std::endl;
    return typename ImageType::Pointer();
  }
};

template<typename PixelType>
struct VoiLoader<PixelType, 3> {
  typedef itk::Image<PixelType, 3> ImageType;

  static typename ImageType::Pointer Load(const char *p_cFileName, const itk::Size<3> &clSize) {
    // VOI
    nih::VoiFile clVoi;
    if (!clVoi.LoadFromFile(p_cFileName))
      return typename ImageType::Pointer();

    typename ImageType::Pointer p_clMask = ImageType::New();
    p_clMask->SetRegions(clSize);
    p_clMask->Allocate(true);

    if (!nih::ConvertVoiToMask3D<PixelType>(p_clMask, clVoi, false))
      return typename ImageType::Pointer();

    return p_clMask;
  }
};

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadImage(const char *p_cFileName);

template<typename PixelType, unsigned int Dimension>
bool SaveImage(typename itk::Image<PixelType, Dimension>::Pointer, const char *p_cFileName);

class CodeBook {
public:
  CodeBook() {
    m_strImageExt = ".mhd";
    m_strForegroundExt = ".voi";
  }

  virtual ~CodeBook() { }

  virtual bool Train(const std::string &strModelFile, const nih::ParameterContainer &clParams) = 0;
  virtual bool LoadModel(const std::string &strModelFile) = 0;
  virtual bool Convert(const std::string &strInputImage, const std::string &strOutputImage) const = 0;
  virtual bool Invert(const std::string &strInputImage, const std::string &strOutputImage) const = 0;
  virtual bool Difference(const std::string &strInputImage, const std::string &strOutputImage) const = 0;
  virtual bool Plot() const = 0;

  virtual void SetImageRoot(const std::string &strImageRoot) {
    m_strImageRoot = strImageRoot;
  }

  virtual void SetImageExt(const std::string &strImageExt) {
    m_strImageExt = strImageExt;
  }

  virtual void SetForegroundRoot(const std::string &strForegroundRoot) {
    m_strForegroundRoot = strForegroundRoot;
  }

  virtual void SetForegroundExt(const std::string &strForegroundExt) {
    m_strForegroundExt = strForegroundExt;
  }

  virtual void SetListFile(const std::string &strListFile) {
    m_strListFile = strListFile;
  }

  const std::string & GetImageRoot() const {
    return m_strImageRoot;
  }

  const std::string & GetImageExt() const {
    return m_strImageExt;
  }

  const std::string & GetForegroundRoot() const {
    return m_strForegroundRoot;
  }

  const std::string & GetForegroundExt() const {
    return m_strForegroundExt;
  }

  const std::string & GetListFile() const {
    return m_strListFile;
  }

  bool LoadList(std::vector<std::string> &vList) const {
    vList.clear();

    std::ifstream listStream(GetListFile().c_str());
    if (!listStream)
      return false;

    std::string strLine;
    while (std::getline(listStream, strLine)) {
      const size_t p = strLine.find('\r');
      if (p != std::string::npos)
        strLine.erase(p);

      vList.push_back(strLine);
    }

    return vList.size() > 0;
  }

private:
  std::string m_strImageRoot, m_strForegroundRoot, m_strListFile, m_strImageExt, m_strForegroundExt;
};

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
class CodeBookTemplate : public CodeBook {
public:
  typedef CodeBook SuperClass;
  typedef CodeBookTree<InputPixelType, OutputPixelType> CodeBookType;
  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
  typedef itk::Image<unsigned char, Dimension> ForegroundMaskType;

  static bool Matches(const std::string &strInputPixelType, const std::string &strOutputPixelType, unsigned int uiDimension) {
    return uiDimension == Dimension && strInputPixelType == GetPixelType<InputPixelType>() && strOutputPixelType == GetPixelType<OutputPixelType>();
  }

  virtual ~CodeBookTemplate() { }

  virtual typename ForegroundMaskType::Pointer LoadForegroundMask(const char *p_cFileName, const itk::Size<Dimension> &clSize) const;

  virtual typename InputImageType::Pointer LoadInputImage(const char *p_cFileName) const {
    return LoadImage<InputPixelType, Dimension>(p_cFileName);
  }

  virtual typename OutputImageType::Pointer LoadOutputImage(const char *p_cFileName) const {
    return LoadImage<OutputPixelType, Dimension>(p_cFileName);
  }

  virtual bool SaveOutputImage(typename OutputImageType::Pointer p_clImage, const char *p_cFileName) const {
    return SaveImage<OutputPixelType, Dimension>(p_clImage, p_cFileName);
  }

  virtual bool SaveInputImage(typename InputImageType::Pointer p_clImage, const char *p_cFileName) const {
    return SaveImage<InputPixelType, Dimension>(p_clImage, p_cFileName);
  }

  virtual bool LoadModel(const std::string &strModelFile) {
    return m_clCodeBook.LoadFromFile(strModelFile.c_str());
  }

  virtual bool ExtractPixels(const std::string &strName, std::vector<InputPixelType> &vPixels) const;

  virtual bool Train(const std::string &strModelFile, const nih::ParameterContainer &clParams);

  virtual bool Convert(const std::string &strInputImage, const std::string &strOutputImage) const;
  virtual bool Invert(const std::string &strInputImage, const std::string &strOutputImage) const;
  virtual bool Difference(const std::string &strInputImage, const std::string &strOutputImage) const;
  virtual bool Plot() const;

private:
  CodeBookType m_clCodeBook;
};

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
class CodeBookTemplateForDicom : public CodeBookTemplate<InputPixelType, OutputPixelType, Dimension> {
public:
  typedef CodeBookTemplate<InputPixelType, OutputPixelType, Dimension> SuperClass;
  typedef typename SuperClass::InputImageType InputImageType;
  typedef typename SuperClass::OutputImageType OutputImageType;
  typedef typename SuperClass::ForegroundMaskType ForegroundMaskType;

  using SuperClass::SetImageExt;
  using SuperClass::GetImageExt;

  static bool Matches(const std::string &strInputPixelType, const std::string &strOutputPixelType, unsigned int uiDimension) {
    return uiDimension == Dimension && strInputPixelType == (std::string("dicom") + GetPixelType<InputPixelType>()) && strOutputPixelType == GetPixelType<OutputPixelType>();
  }

  CodeBookTemplateForDicom() {
    SetImageExt("");
  }

  virtual ~CodeBookTemplateForDicom() { }

  virtual typename InputImageType::Pointer LoadInputImage(const char *p_cFileName) const {
    if (Dimension != 2 && Dimension != 3)
      return typename InputImageType::Pointer();

    if (Dimension == 2)
      return nih::LoadDicomImage<InputPixelType, Dimension>(p_cFileName);

    if (nih::IsFolder(p_cFileName) || nih::FileExists(p_cFileName))
      return nih::LoadDicomImage<InputPixelType, Dimension>(p_cFileName);

    // Try Prostate convention
    std::string strPattern = p_cFileName;

    if (strPattern.find("%d") == std::string::npos) {
      strPattern += ".%d.dcm";
      strPattern += GetImageExt();
    }

    return nih::LoadDicomImageByIndex<InputPixelType, Dimension>(strPattern.c_str());
  }

  virtual bool SaveInputImage(typename InputImageType::Pointer p_clImage, const char *p_cFileName) const {
    std::cerr << "Warning: Saving DICOM sequences is not supported. Saving as a non-DICOM image instead." << std::endl;
    return SaveImage<InputPixelType, Dimension>(p_clImage, p_cFileName);
  }
};

class CodeBookFactory {
public:
  class Node {
  public:
    virtual ~Node() { }
    virtual bool Matches(const std::string &strInputPixelType, const std::string &strOutputPixelType, unsigned int uiDimension) = 0;
    virtual std::shared_ptr<CodeBook> Create() = 0;
  };

  template<typename CodeBookType>
  class NodeTemplate : public Node {
  public:
    virtual ~NodeTemplate() { }
    virtual bool Matches(const std::string &strInputPixelType, const std::string &strOutputPixelType, unsigned int uiDimension) {
      return CodeBookType::Matches(strInputPixelType, strOutputPixelType, uiDimension);
    }
    virtual std::shared_ptr<CodeBook> Create() {
      return std::make_shared<CodeBookType>();
    }
  };

  static CodeBookFactory & GetInstance() {
    static CodeBookFactory clFactory;
    return clFactory;
  }

  template<typename CodeBookType>
  void Register() {
    m_vNodes.push_back(std::make_shared<NodeTemplate<CodeBookType> >());
  }

  std::shared_ptr<CodeBook> Create(const std::string &strInputPixelType, const std::string &strOutputPixelType, unsigned int uiDimension) const {
    for (size_t i = 0; i < m_vNodes.size(); ++i) {
      if (m_vNodes[i]->Matches(strInputPixelType, strOutputPixelType, uiDimension))
        return m_vNodes[i]->Create();
    }

    return std::shared_ptr<CodeBook>();
  }

private:
  std::vector<std::shared_ptr<Node> > m_vNodes;

  CodeBookFactory() {
    // 2D
    Register<CodeBookTemplate<ushort, uchar, 2> >();
    Register<CodeBookTemplate<short, uchar, 2> >();
    Register<CodeBookTemplate<float, uchar, 2> >();
    Register<CodeBookTemplate<double, uchar, 2> >();
    Register<CodeBookTemplate<float, ushort, 2> >();
    Register<CodeBookTemplate<double, ushort, 2> >();
    Register<CodeBookTemplateForDicom<ushort, uchar, 2> >();
    Register<CodeBookTemplateForDicom<short, uchar, 2> >();

    // 3D
    Register<CodeBookTemplate<ushort, uchar, 3> >();
    Register<CodeBookTemplate<short, uchar, 3> >();
    Register<CodeBookTemplate<float, uchar, 3> >();
    Register<CodeBookTemplate<double, uchar, 3> >();
    Register<CodeBookTemplate<float, ushort, 3> >();
    Register<CodeBookTemplate<double, ushort, 3> >();
    Register<CodeBookTemplateForDicom<ushort, uchar, 3> >();
    Register<CodeBookTemplateForDicom<short, uchar, 3> >();
  }
};

int main(int argc, char **argv) {
  const char * const p_cArg0 = argv[0];

  if (argc < 2)
    Usage(p_cArg0);

  std::string strMode = argv[1];
  unsigned int uiDimension = 3;
  double dMinSampleSize = 1000.0;
  std::string strImageExt;
  std::string strCodeBookFile = "codebook.clf";
  std::string strInputPixelType = "short";
  std::string strOutputPixelType = "uchar";
  std::string strListFile;
  std::string strMaskExt;
  std::string strMaskRoot;
  std::string strImageRoot;

  optind = 2;

  int c = 0;
  while ((c = getopt(argc, argv, "d:e:f:hi:k:l:m:o:r:s:")) != -1) {
    switch (c) {
    case 'd':
      {
        char *p = NULL;
        uiDimension = strtoul(optarg, &p, 10);

        if (*p != '\0')
          Usage(p_cArg0);
      }
      break;
    case 'e':
      strImageExt = optarg;
      break;
    case 'f':
      strCodeBookFile = optarg;
      break;
    case 'h':
      Usage(p_cArg0);
      break;
    case 'i':
      strInputPixelType = optarg;
      break;
    case 'k':
      strMaskRoot = optarg;
      break;
    case 'l':
      strListFile = optarg;
      break;
    case 'm': 
      strMaskExt = optarg;
      break;
    case 'o':
      strOutputPixelType = optarg;
      break;
    case 'r':
      strImageRoot = optarg;
      break;
    case 's':
      {
        char *p = NULL;
        dMinSampleSize = strtod(optarg, &p);
        if (*p != '\0')
          Usage(p_cArg0);
      }
      break;
    case '?':
    default:
      Usage(p_cArg0);
    }
  }

  argc -= optind;
  argv += optind;

  if (strMode != "train" && strMode != "convert" && strMode != "invert" && strMode != "plot" && strMode != "difference") {
    std::cerr << "Error: Invalid mode '" << strMode << "'." << std::endl;
    Usage(p_cArg0);
  }

  if (strInputPixelType.empty() || strOutputPixelType.empty()) {
    std::cerr << "Error: Input/Output pixel types not specified." << std::endl;
    return -1;
  }

  CodeBookFactory &clFactory = CodeBookFactory::GetInstance();

  std::shared_ptr<CodeBook> p_clCodeBook = clFactory.Create(strInputPixelType, strOutputPixelType, uiDimension);
  if (!p_clCodeBook) {
    std::cerr << "Error: Could not create CodeBook." << std::endl;
    return -1;
  }

  if (strMode == "train") {
    if (strListFile.empty()) {
      std::cerr << "Error: No list file specified." << std::endl;
      return -1;
    }

    if (strImageRoot.empty()) {
      std::cerr << "Error: No image root specified." << std::endl;
      return -1;
    }

    p_clCodeBook->SetListFile(strListFile);
    p_clCodeBook->SetImageRoot(strImageRoot);

    if (strMaskRoot.size() > 0)
      p_clCodeBook->SetForegroundRoot(strMaskRoot);

    if (strImageExt.size() > 0)
      p_clCodeBook->SetImageExt(strImageExt);

    if (strMaskExt.size() > 0)
      p_clCodeBook->SetForegroundExt(strMaskExt);

    nih::ParameterMap clParams;
    clParams.SetValue("minSampleSize", dMinSampleSize);

    if (!p_clCodeBook->Train(strCodeBookFile, clParams)) {
      std::cerr << "Error: Training failed." << std::endl;
      return -1;
    }
  }
  else if (strMode == "convert") {
    if (argc != 2)
      Usage(p_cArg0);

    const std::string strInputImage = argv[0];
    const std::string strOutputImage = argv[1];

    if (!p_clCodeBook->LoadModel(strCodeBookFile)) {
      std::cerr << "Error: Failed to load code book '" << strCodeBookFile << "'." << std::endl;
      return -1;
    }

    if (!p_clCodeBook->Convert(strInputImage, strOutputImage)) {
      std::cerr << "Error: Convert failed." << std::endl;
      return -1;
    }
  }
  else if (strMode == "invert") {
    if (argc != 2)
      Usage(p_cArg0);

    const std::string strInputImage = argv[0];
    const std::string strOutputImage = argv[1];

    if (!p_clCodeBook->LoadModel(strCodeBookFile)) {
      std::cerr << "Error: Failed to load code book '" << strCodeBookFile << "'." << std::endl;
      return -1;
    }

    if (!p_clCodeBook->Invert(strInputImage, strOutputImage)) {
      std::cerr << "Error: Convert failed." << std::endl;
      return -1;
    }
  }
  else if (strMode == "difference") {
    if (argc != 2)
      Usage(p_cArg0);

    const std::string strInputImage = argv[0];
    const std::string strOutputImage = argv[1];

    if (!p_clCodeBook->LoadModel(strCodeBookFile)) {
      std::cerr << "Error: Failed to load code book '" << strCodeBookFile << "'." << std::endl;
      return -1;
    }

    if (!p_clCodeBook->Difference(strInputImage, strOutputImage)) {
      std::cerr << "Error: Convert failed." << std::endl;
      return -1;
    }
  }
  else if (strMode == "plot") {
    if (!p_clCodeBook->LoadModel(strCodeBookFile)) {
      std::cerr << "Error: Failed to load code book '" << strCodeBookFile << "'." << std::endl;
      return -1;
    }

    if (!p_clCodeBook->Plot()) {
      std::cerr << "Error: Failed to plot pixel values." << std::endl;
      return -1;
    }
  }

  return 0;
}

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadImage(const char *p_cFileName) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer p_clReader = ReaderType::New();

  p_clReader->SetFileName(p_cFileName);

  try {
    p_clReader->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  return p_clReader->GetOutput();
}

template<typename PixelType, unsigned int Dimension>
bool SaveImage(typename itk::Image<PixelType, Dimension>::Pointer p_clImage, const char *p_cFileName) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  if (!p_clImage)
    return false;

  typename WriterType::Pointer p_clWriter = WriterType::New();

  p_clWriter->SetFileName(p_cFileName);
  p_clWriter->SetInput(p_clImage);

  try {
    p_clWriter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return false;
  }

  return true;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
typename CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::ForegroundMaskType::Pointer CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::LoadForegroundMask(const char *p_cFileName, const itk::Size<Dimension> &clSize) const {
  typename ForegroundMaskType::Pointer p_clMask;

  if (Dimension == 3 && strstr(p_cFileName, ".voi") != NULL)
    return VoiLoader<unsigned char, Dimension>::Load(p_cFileName, clSize);

  p_clMask = LoadImage<unsigned char, Dimension>(p_cFileName);

  if (!p_clMask || p_clMask->GetBufferedRegion().GetSize() != clSize)
    return typename ForegroundMaskType::Pointer();

  return p_clMask;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::ExtractPixels(const std::string &strName, std::vector<InputPixelType> &vPixels) const {
  const std::string strImageFileName = GetImageRoot() + '/' + strName + GetImageExt();

  typename InputImageType::Pointer p_clImage = LoadInputImage(strImageFileName.c_str());
  if (!p_clImage) {
    std::cerr << "Error: Could not load image '" << strImageFileName << "'." << std::endl;
    return false;
  }

  const InputPixelType * const p_inBuffer = p_clImage->GetPixelContainer()->GetBufferPointer();
  const size_t length = p_clImage->GetPixelContainer()->Size();
  const unsigned char * p_ucMaskBuffer = NULL;

  typename ForegroundMaskType::Pointer p_clMask;

  if (GetForegroundRoot().size() > 0) {
    const std::string strForegroundFileName = GetForegroundRoot() + '/' + strName + GetForegroundExt();

    p_clMask = LoadForegroundMask(strForegroundFileName.c_str(), p_clImage->GetBufferedRegion().GetSize());
    if (!p_clMask) {
      std::cerr << "Error: Could not load foreground '" << strForegroundFileName << "'." << std::endl;
      return false;
    }

    p_ucMaskBuffer = p_clMask->GetPixelContainer()->GetBufferPointer();
  }

  for (size_t i = 0; i < length; ++i) {
    const bool bInForeground = p_ucMaskBuffer != NULL ? (p_ucMaskBuffer[i] != 0) : true;

    if (bInForeground)
      vPixels.push_back(p_inBuffer[i]);
  }

  return true;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::Train(const std::string &strModelFile, const nih::ParameterContainer &clParams) {
  double dMinSampleSize = clParams.GetValue<double>("minSampleSize", 1000.0);

  if (dMinSampleSize < 0.0)
    return false;

  m_clCodeBook.Clear();


  std::vector<std::string> vList;
  if (!LoadList(vList)) {
    std::cerr << "Error: Could not load training list." << std::endl;
    return false;
  }

  std::vector<InputPixelType> vPixels;

  for (size_t i = 0; i < vList.size(); ++i) {
    if (!ExtractPixels(vList[i], vPixels)) {
      std::cerr << "Error: Failed to extract pixels." << std::endl;
      return false;
    }
  }

  if (dMinSampleSize > 0.0 && dMinSampleSize < 1.0)
    m_clCodeBook.SetMinSampleSize((unsigned int)(dMinSampleSize * vPixels.size() + 0.5f));
  else
    m_clCodeBook.SetMinSampleSize((unsigned int)dMinSampleSize);

  std::cout << "Info: Training sample size = " << vPixels.size() << std::endl;

  //m_clCodeBook.SetMinSampleSize(vPixels.size() / m_clCodeBook.GetMaxNumberOfLeaves());

  if (!m_clCodeBook.Train(vPixels)) {
    std::cerr << "Error: Failed to train codebook." << std::endl;
    return false;
  }

  std::cout << "Info: Number of codes used = " << m_clCodeBook.CountLeaves() << std::endl;

  return m_clCodeBook.SaveToFile(strModelFile.c_str());
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::Convert(const std::string &strInputImage, const std::string &strOutputImage) const {
  if (m_clCodeBook.Empty()) {
    std::cerr << "Error: No CodeBook loaded." << std::endl;
    return false;
  }

  typename InputImageType::Pointer p_clInputImage = LoadInputImage(strInputImage.c_str());

  if (!p_clInputImage) {
    std::cerr << "Error: Could not load image '" << strInputImage << "'."  << std::endl;
    return false;
  }

  typename OutputImageType::Pointer p_clOutputImage = OutputImageType::New();

  p_clOutputImage->SetRegions(p_clInputImage->GetBufferedRegion().GetSize());
  p_clOutputImage->SetSpacing(p_clInputImage->GetSpacing());
  p_clOutputImage->SetOrigin(p_clInputImage->GetOrigin());
  p_clOutputImage->SetDirection(p_clInputImage->GetDirection());
  p_clOutputImage->SetMetaDataDictionary(p_clInputImage->GetMetaDataDictionary());

  p_clOutputImage->Allocate();

  const InputPixelType * const p_inBuffer = p_clInputImage->GetPixelContainer()->GetBufferPointer();
  const size_t length = p_clInputImage->GetPixelContainer()->Size();
  OutputPixelType * const p_outBuffer = p_clOutputImage->GetPixelContainer()->GetBufferPointer();

  for (size_t i = 0; i < length; ++i) {
    if (!m_clCodeBook.Convert(p_inBuffer[i], p_outBuffer[i])) {
      std::cerr << "Error: Conversion failed." << std::endl;
      return false;
    }
  }

  if (!SaveOutputImage(p_clOutputImage, strOutputImage.c_str())) {
    std::cerr << "Error: Could not save '" << strOutputImage << "'." << std::endl;
    return false;
  }

  return true;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::Invert(const std::string &strInputImage, const std::string &strOutputImage) const {
  if (m_clCodeBook.Empty()) {
    std::cerr << "Error: No CodeBook loaded." << std::endl;
    return false;
  }

  typename OutputImageType::Pointer p_clInputImage = LoadOutputImage(strInputImage.c_str());

  if (!p_clInputImage) {
    std::cerr << "Error: Could not load image '" << strInputImage << "'."  << std::endl;
    return false;
  }

  typename InputImageType::Pointer p_clOutputImage = InputImageType::New();

  p_clOutputImage->SetRegions(p_clInputImage->GetBufferedRegion().GetSize());
  p_clOutputImage->SetSpacing(p_clInputImage->GetSpacing());
  p_clOutputImage->SetOrigin(p_clInputImage->GetOrigin());
  p_clOutputImage->SetDirection(p_clInputImage->GetDirection());
  p_clOutputImage->SetMetaDataDictionary(p_clInputImage->GetMetaDataDictionary());

  p_clOutputImage->Allocate();

  const OutputPixelType * const p_inBuffer = p_clInputImage->GetPixelContainer()->GetBufferPointer();
  const size_t length = p_clInputImage->GetPixelContainer()->Size();
  InputPixelType * const p_outBuffer = p_clOutputImage->GetPixelContainer()->GetBufferPointer();

  for (size_t i = 0; i < length; ++i) {
    if (!m_clCodeBook.Invert(p_inBuffer[i], p_outBuffer[i])) {
      std::cerr << "Error: Conversion failed." << std::endl;
      return false;
    }
  }

  if (!SaveInputImage(p_clOutputImage, strOutputImage.c_str())) {
    std::cerr << "Error: Could not save '" << strOutputImage << "'." << std::endl;
    return false;
  }

  return true;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::Difference(const std::string &strInputImage, const std::string &strOutputImage) const {
  if (m_clCodeBook.Empty()) {
    std::cerr << "Error: No CodeBook loaded." << std::endl;
    return false;
  }

  typename InputImageType::Pointer p_clInputImage = LoadInputImage(strInputImage.c_str());

  if (!p_clInputImage) {
    std::cerr << "Error: Could not load image '" << strInputImage << "'."  << std::endl;
    return false;
  }

  typename InputImageType::Pointer p_clOutputImage = InputImageType::New();

  p_clOutputImage->SetRegions(p_clInputImage->GetBufferedRegion().GetSize());
  p_clOutputImage->SetSpacing(p_clInputImage->GetSpacing());
  p_clOutputImage->SetOrigin(p_clInputImage->GetOrigin());
  p_clOutputImage->SetDirection(p_clInputImage->GetDirection());
  p_clOutputImage->SetMetaDataDictionary(p_clInputImage->GetMetaDataDictionary());

  p_clOutputImage->Allocate();

  const InputPixelType * const p_inBuffer = p_clInputImage->GetPixelContainer()->GetBufferPointer();
  const size_t length = p_clInputImage->GetPixelContainer()->Size();
  InputPixelType * const p_outBuffer = p_clOutputImage->GetPixelContainer()->GetBufferPointer();

  double dRMSE= 0.0;

  for (size_t i = 0; i < length; ++i) {
    OutputPixelType outPixel = OutputPixelType();
    if (!m_clCodeBook.Convert(p_inBuffer[i], outPixel) || !m_clCodeBook.Invert(outPixel, p_outBuffer[i])) {
      std::cerr << "Error: Conversion failed." << std::endl;
      return false;
    }

    p_outBuffer[i] -= p_inBuffer[i];

    dRMSE += std::pow((double)p_outBuffer[i], 2);
  }

  dRMSE = length > 0 ? std::sqrt(dRMSE/length) : -1.0;

  std::cout << "RMSE: " << dRMSE << std::endl;

  if (!SaveInputImage(p_clOutputImage, strOutputImage.c_str())) {
    std::cerr << "Error: Could not save '" << strOutputImage << "'." << std::endl;
    return false;
  }

  return true;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
bool CodeBookTemplate<InputPixelType, OutputPixelType, Dimension>::Plot() const {
  if (m_clCodeBook.Empty())
    return false;

  const InputPixelType xBegin = m_clCodeBook.MinPixel();
  const InputPixelType xEnd = m_clCodeBook.MaxPixel();
  const InputPixelType xStep = StepValue(xBegin, xEnd);

  for (InputPixelType x = xBegin; x <= xEnd; x += xStep) {
    const InputPixelType inPixel = x;
    OutputPixelType outPixel = OutputPixelType();
    if (!m_clCodeBook.Convert(inPixel, outPixel)) {
      std::cerr << "Error: Failed to convert pixel value " << x << std::endl;
      return false;
    }

    std::cout << x << ' ' << (int64_t)outPixel << std::endl;

    if (x == std::numeric_limits<InputPixelType>::max()) // Prevent overflow!
      break;
  }

  return true;
}
