#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_set>
#include "itkImage.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMImageIOFactory.h"
#include "itkMetaImageIO.h"
#include "itkMetaImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkNiftiImageIO.h"
#include "itkNiftiImageIOFactory.h"
#include "itkNrrdImageIO.h"
#include "itkNrrdImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

#ifdef _WIN32
#undef CopyFile
#undef MoveFile
#endif // _WIN32

#define MAX_NUM_THREADS 32

namespace nih {

template<typename DataType>
struct NoDeleter {
  void operator()(DataType *) const { }
};

template<typename PixelType>
struct NonMaximumSuppressionCompare {
  bool operator()(const PixelType *a, const PixelType *b) const {
    return *a > *b;
  }
};

template<typename RealType>
RealType Round(const RealType &x) {
  RealType intPart = RealType();
  std::modf(std::abs(x) + RealType(0.5), &intPart);
  return x < 0 ? -intPart : intPart;  
}

template<typename IteratorType>
void PartialShuffle(IteratorType begin, IteratorType stop, IteratorType end) {
  size_t length = std::distance(begin,end);

  while (begin != stop && length > 1) {
#ifdef _WIN32
    const uint64_t r = (uint64_t)rand() | ((uint64_t)rand() << 15) | ((uint64_t)rand() << 30) | ((uint64_t)rand() << 45);
#else // !_WIN32
    const uint64_t r = (uint64_t)rand() | ((uint64_t)rand() << 31);
#endif
    IteratorType itr(begin);
    std::advance(itr, (r % length));
    std::iter_swap(begin, itr);
    --length;
    ++begin;
  }
}

template<typename PixelType, unsigned int Dimension>
class ImageAccessor {
public:
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::Point<float, Dimension> WorldPointType;
  typedef itk::Point<int, Dimension> PixelPointType;

  typename ImageType::Pointer p_clImage;

  ImageAccessor() { }
  ImageAccessor(typename ImageType::Pointer p_clImage_)
  : p_clImage(p_clImage_) { }

  bool Valid(const PixelPointType &clPoint) const {
    if (!p_clImage)
      return false;

    typedef typename ImageType::SizeType SizeType;
    const SizeType &clSize = p_clImage->GetBufferedRegion().GetSize();

    for (unsigned int i = 0; i < Dimension; ++i) {
      if (clPoint[i] < 0 || clPoint[i] >= clSize[i])
        return false;
    }

    return true;
  }

  bool Valid(const WorldPointType &clPointMM) const {
    return Valid(ToPixel(clPointMM));
  }

  PixelPointType ToPixel(const WorldPointType &clPointMM) const {
    if (!p_clImage)
      return PixelPointType(0);

    typedef typename ImageType::SpacingType SpacingType;
    typedef typename ImageType::PointType OriginPointType;

    const OriginPointType &clOrigin = p_clImage->GetOrigin();
    const SpacingType &clSpacing = p_clImage->GetSpacing();

    PixelPointType clPoint;
    for (unsigned int i = 0; i < Dimension; ++i)
      clPoint[i] = (int)Round((clPointMM[i] - clOrigin[i])/clSpacing[i]);

    return clPoint;
  }

  WorldPointType ToWorld(const PixelPointType &clPoint) const {
    if (!p_clImage)
      return WorldPointType(0.0f);

    typedef typename ImageType::SpacingType SpacingType;
    typedef typename ImageType::PointType OriginPointType;

    const OriginPointType &clOrigin = p_clImage->GetOrigin();
    const SpacingType &clSpacing = p_clImage->GetSpacing();

    WorldPointType clPointMM;
    for (unsigned int i = 0; i < Dimension; ++i)
      clPointMM[i] = (float)(clOrigin[i] + clSpacing[i]*clPoint[i]);

    return clPointMM;
  }

  PixelType GetPixel(const PixelPointType &clPoint) const {
    if (!p_clImage)
      return PixelType();

    typedef typename ImageType::IndexType IndexType;
    IndexType clIndex;
    for (unsigned int i = 0; i < Dimension; ++i)
      clIndex[i] = clPoint[i];

    return p_clImage->GetPixel(clIndex);
  }

  PixelType GetPixel(const WorldPointType &clPointMM) const {
    return GetPixel(ToPixel(clPointMM));
  }

  void SetPixel(const PixelPointType &clPoint, const PixelType &value) const {
    if (!p_clImage)
      return;

    typedef typename ImageType::IndexType IndexType;
    IndexType clIndex;
    for (unsigned int i = 0; i < Dimension; ++i)
      clIndex[i] = clPoint[i];

    p_clImage->SetPixel(clIndex, value);
  }

  void SetPixel(const WorldPointType &clPointMM, const PixelType &value) const {
    SetPixel(ToPixel(clPointMM), value);
  }
};

int64_t Rand64();

bool IsFolder(const char *p_cPath);
bool MoveFile(const char *p_cSource, const char *p_cDest);
bool CopyFile(const char *p_cSource, const char *p_cDest);
bool MakeFolder(const char *p_cPath);
bool DeleteFolder(const char *p_cPath);
bool FileExists(const char *p_cPath);
void FindDicomFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder);
void FindFiles(std::vector<std::string> &vFiles, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive = false);
void FindFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive = false);

std::string GetCurrentTimeString();

std::string DirName(const std::string &strFilePath);
std::string BaseName(const std::string &strFilePath);
std::string StripExtension(const std::string &strFilePath);
inline std::string ReplaceExtension(const std::string &strFilePath, const std::string &strExt) {
  return StripExtension(strFilePath) + strExt;
}

void Trim(std::string &strString);

void SetLogFile(const char *p_cLogFile);
void ClearLogFile();
void Log(const char *p_cFormat, ...);

void RegisterImageFormats();

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
typename itk::Image<OutputPixelType, Dimension>::Pointer GaussianSmoothing(typename itk::Image<InputPixelType, Dimension>::Pointer p_clImage, const float *p_fVar);

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer CopyImage(typename itk::Image<PixelType, Dimension>::Pointer p_clImage);

template<typename PixelType>
typename itk::Image<PixelType,3>::Pointer NonMaximumSuppression3D(typename itk::Image<PixelType,3>::Pointer p_clImage, unsigned int uiWindowWidth, unsigned int uiWindowHeight, unsigned int uiWindowDepth);

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,3>::Pointer LabelConnectedComponents3D(typename itk::Image<InputPixelType,3>::Pointer p_clMask, unsigned int &uiNumComponents);

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,2>::Pointer LabelConnectedComponents2D(typename itk::Image<InputPixelType,2>::Pointer p_clMask, unsigned int &uiNumComponents);

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
std::pair<typename itk::Image<OutputPixelType,Dimension>::Pointer, typename itk::Image<InputPixelType,Dimension>::Pointer> ComputeDistanceMapWithLabels(typename itk::Image<InputPixelType,Dimension>::Pointer p_clMask);

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
typename itk::Image<OutputPixelType,Dimension>::Pointer ComputeDistanceMap(typename itk::Image<InputPixelType,Dimension>::Pointer p_clMask) {
  return ComputeDistanceMapWithLabels<InputPixelType,OutputPixelType,Dimension>(p_clMask).first;
}

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,2>::Pointer ComputeDistanceMap2D(typename itk::Image<InputPixelType,2>::Pointer p_clMask) {
  return ComputeDistanceMap<InputPixelType,OutputPixelType,2>(p_clMask);
}

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,3>::Pointer ComputeDistanceMap3D(typename itk::Image<InputPixelType,3>::Pointer p_clMask) {
  return ComputeDistanceMap<InputPixelType,OutputPixelType,3>(p_clMask);
}

template<typename PixelType>
bool GetSlice(typename itk::Image<PixelType, 2>::Pointer p_clSlice, typename itk::Image<PixelType, 3>::Pointer p_clVolume, unsigned int z);

template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer GetSlice(typename itk::Image<PixelType, 3>::Pointer p_clVolume, unsigned int z) {
  typedef itk::Image<PixelType, 2> ImageType2D;

  typename ImageType2D::Pointer p_clSlice = ImageType2D::New();

  if (!GetSlice<PixelType>(p_clSlice, p_clVolume, z))
    return typename ImageType2D::Pointer();

  return p_clSlice;
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer LoadVolume(const char *p_cFileName);

template<typename PixelType>
bool SaveVolume(typename itk::Image<PixelType, 3>::Pointer p_clVolume, const char *p_cFileName, bool bCompress = true);

template<typename PixelType>
bool SaveSlice(typename itk::Image<PixelType, 2>::Pointer p_clImage, const char *p_cFileName, bool bCompress = true);

template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer ResampleSlice(typename itk::Image<PixelType, 2>::Pointer p_clSlice, int iNewWidth, int iNewHeight, bool bNearest = false);

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer ResampleSlices(typename itk::Image<PixelType, 3>::Pointer p_clVolume, int iNewWidth, int iNewHeight, bool bNearest = false);

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer ResampleVolume(typename itk::Image<PixelType, 3>::Pointer p_clVolume, int iNewWidth, int iNewHeight, int iNewDepth, bool bNearest = false);

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadDicomImage(const char *p_cPath, const char *p_cSeries = NULL);

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadDicomImageByIndex(const char *p_cPattern, int iBeginIndex = 0, int iEndIndex = -1);

template<typename PixelType>
typename itk::Image<PixelType, 4>::Pointer LoadDicomVideo(const char *p_cPath, const char *p_cSeries = NULL);

inline itk::Image<short, 3>::Pointer LoadDicomVolume(const char *p_cPath, const char *p_cSeries = NULL) {
  return LoadDicomImage<short, 3>(p_cPath, p_cSeries);
}

inline itk::Image<short, 3>::Pointer LoadDicomVolumeByIndex(const char *p_cPattern, int iBeginIndex = 0, int iEndIndex = -1) {
  return LoadDicomImageByIndex<short, 3>(p_cPattern, iBeginIndex, iEndIndex);
}

inline itk::Image<short, 2>::Pointer LoadDicomSlice(const char *p_cPath) {
  return LoadDicomImage<short, 2>(p_cPath);
}

template<typename PixelType>
bool SaveDicomSliceT(typename itk::Image<PixelType, 2>::Pointer p_clSlice, const char *p_cFileName, bool bCompress = false);

inline bool SaveDicomSlice(itk::Image<short, 2>::Pointer p_clSlice, const char *p_cFileName, bool bCompress = false) {
  return SaveDicomSliceT<short>(p_clSlice, p_cFileName, bCompress);
}

bool GetBValueFileNames(std::map<double, std::vector<std::string> > &mBValueMap, const char * p_cPath);
itk::Image<short, 3>::Pointer LoadB2000Volume(const char *p_cPath, bool bNormalize = true);

template<unsigned int Dimension>
bool NormalizeB2000(typename itk::Image<short,Dimension>::Pointer p_clB2000Image, typename itk::Image<short,Dimension>::Pointer p_clB0Image);

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
typename itk::Image<OutputPixelType, Dimension>::Pointer GaussianSmoothing(typename itk::Image<InputPixelType, Dimension>::Pointer p_clImage, const float *p_fVar) {
  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
  typedef itk::DiscreteGaussianImageFilter<InputImageType, OutputImageType> FilterType;

  if (!p_clImage)
    return typename OutputImageType::Pointer();

  typename FilterType::Pointer p_clFilter = FilterType::New();

  p_clFilter->SetVariance(p_fVar);
  p_clFilter->SetUseImageSpacing(true);
  p_clFilter->SetInput(p_clImage);

  try {
    p_clFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return typename OutputImageType::Pointer();
  }

  return p_clFilter->GetOutput();
}

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer CopyImage(typename itk::Image<PixelType, Dimension>::Pointer p_clImage) {
  typedef itk::Image<PixelType, Dimension> ImageType;

  typename ImageType::Pointer p_clOutputImage = ImageType::New();

  if (!p_clOutputImage)
    return typename ImageType::Pointer();

  p_clOutputImage->SetRegions(p_clImage->GetBufferedRegion());
  p_clOutputImage->SetSpacing(p_clImage->GetSpacing());
  p_clOutputImage->SetOrigin(p_clImage->GetOrigin());
  p_clOutputImage->SetMetaDataDictionary(p_clImage->GetMetaDataDictionary());

  p_clOutputImage->Allocate(false);

  const size_t bufferSize = p_clImage->GetPixelContainer()->Size();

  std::memcpy(p_clOutputImage->GetPixelContainer()->GetImportPointer(), p_clImage->GetPixelContainer()->GetImportPointer(), bufferSize*sizeof(PixelType));

  return p_clOutputImage;
}

template<typename PixelType>
typename itk::Image<PixelType,3>::Pointer NonMaximumSuppression3D(typename itk::Image<PixelType,3>::Pointer p_clImage, unsigned int uiWindowWidth, unsigned int uiWindowHeight, unsigned int uiWindowDepth) {
  typedef itk::Image<PixelType,3> ImageType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;

  if (!p_clImage)
    return typename ImageType::Pointer();

  typename ImageType::Pointer p_clOutputImage = CopyImage<PixelType, 3>(p_clImage);

  if (!p_clOutputImage)
    return typename ImageType::Pointer();

  const SizeType &clSize = p_clOutputImage->GetBufferedRegion().GetSize();

  const PixelType * const p_buffer = p_clOutputImage->GetPixelContainer()->GetImportPointer();
  const size_t bufferSize = p_clOutputImage->GetPixelContainer()->Size();

  std::vector<const PixelType *> vPointScores;

  for (size_t i = 0; i < bufferSize; ++i) {
    if (p_buffer[i] > 0)
      vPointScores.push_back(p_buffer + i);
  }

  std::sort(vPointScores.begin(), vPointScores.end(), NonMaximumSuppressionCompare<PixelType>());

  // XXX: 0 is assumed to be lowest value
  for (int i = 0; (size_t)i < vPointScores.size(); ++i) {
    const PixelType * const p_pixel = vPointScores[i];

    if (*p_pixel <= 0)
      continue;

    const PixelType score = *p_pixel;

    const IndexType clTmp = p_clOutputImage->ComputeIndex(p_pixel - p_buffer);

    const int x = (int)clTmp[0];
    const int y = (int)clTmp[1];
    const int z = (int)clTmp[2];

    int xBegin = x - (int)uiWindowWidth/2;
    int xEnd = xBegin + (int)uiWindowWidth;

    xBegin = std::max(0, xBegin);
    xEnd = std::min(xEnd, (int)clSize[0]);

    int yBegin = y - (int)uiWindowHeight/2;
    int yEnd = yBegin + (int)uiWindowHeight;

    yBegin = std::max(0, yBegin);
    yEnd = std::min(yEnd, (int)clSize[1]);

    int zBegin = z - (int)uiWindowDepth/2;
    int zEnd = zBegin + (int)uiWindowDepth;

    zBegin = std::max(0, zBegin);
    zEnd = std::min(zEnd, (int)clSize[2]);

    for (int zi = zBegin; zi < zEnd; ++zi) {
      for (int yi = yBegin; yi < yEnd; ++yi) {
        for (int xi = xBegin; xi < xEnd; ++xi) {
          const IndexType clIndex = { xi, yi, zi };

          if (p_clOutputImage->GetPixel(clIndex) <= score)
            p_clOutputImage->SetPixel(clIndex, 0);
        }
      }
    }

    p_clOutputImage->SetPixel(clTmp, score);
  }
  
  return p_clOutputImage;
}

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,3>::Pointer LabelConnectedComponents3D(typename itk::Image<InputPixelType,3>::Pointer p_clMask, unsigned int &uiNumComponents) {
  typedef itk::Image<OutputPixelType, 3> OutputImageType;
  typedef itk::Image<InputPixelType, 3> InputImageType;
  typedef itk::ConnectedComponentImageFilter<InputImageType, OutputImageType> ConnectedComponentFilterType;

  typename ConnectedComponentFilterType::Pointer p_clCCFilter = ConnectedComponentFilterType::New();
  if (!p_clCCFilter)
    return typename OutputImageType::Pointer();

  p_clCCFilter->SetInput(p_clMask);

  try {
    p_clCCFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return typename OutputImageType::Pointer();
  }

  uiNumComponents = (unsigned int)p_clCCFilter->GetObjectCount();
  std::cout << "Info: Number of connected components: " << p_clCCFilter->GetObjectCount() << std::endl;

  typename OutputImageType::Pointer p_clOutputImage = p_clCCFilter->GetOutput();

  return p_clOutputImage;
}

template<typename InputPixelType, typename OutputPixelType>
typename itk::Image<OutputPixelType,2>::Pointer LabelConnectedComponents2D(typename itk::Image<InputPixelType,2>::Pointer p_clMask, unsigned int &uiNumComponents) {
  typedef itk::Image<OutputPixelType, 2> OutputImageType;
  typedef itk::Image<InputPixelType, 2> InputImageType;
  typedef itk::ConnectedComponentImageFilter<InputImageType, OutputImageType> ConnectedComponentFilterType;

  typename ConnectedComponentFilterType::Pointer p_clCCFilter = ConnectedComponentFilterType::New();
  if (!p_clCCFilter)
    return typename OutputImageType::Pointer();

  p_clCCFilter->SetInput(p_clMask);

  try {
    p_clCCFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return typename OutputImageType::Pointer();
  }

  uiNumComponents = (unsigned int)p_clCCFilter->GetObjectCount();
  std::cout << "Info: Number of connected components: " << p_clCCFilter->GetObjectCount() << std::endl;

  typename OutputImageType::Pointer p_clOutputImage = p_clCCFilter->GetOutput();

  return p_clOutputImage;
}

template<typename InputPixelType, typename OutputPixelType, unsigned int Dimension>
std::pair<typename itk::Image<OutputPixelType,Dimension>::Pointer, typename itk::Image<InputPixelType,Dimension>::Pointer> ComputeDistanceMapWithLabels(typename itk::Image<InputPixelType,Dimension>::Pointer p_clMask) {
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef std::pair<typename OutputImageType::Pointer, typename InputImageType::Pointer> PairType;
  typedef itk::SignedDanielssonDistanceMapImageFilter<InputImageType, OutputImageType> DistanceMapFilterType;

  if (!p_clMask)
    return PairType();

  typename DistanceMapFilterType::Pointer p_clDistanceMapFilter = DistanceMapFilterType::New();

  if (!p_clDistanceMapFilter)
    return PairType();

  // ITK: The inside is considered as having negative distances. Outside is treated as having positive distances. To change the convention, use the InsideIsPositive(bool) function.

  p_clDistanceMapFilter->SetInput(p_clMask);
  p_clDistanceMapFilter->SetSquaredDistance(false);
  p_clDistanceMapFilter->SetUseImageSpacing(true);
  p_clDistanceMapFilter->SetInsideIsPositive(false); // Enforce the comment

  //p_clDistanceMapFilter->SetNumberOfThreads(1);

  try {
    p_clDistanceMapFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return PairType();
  }

  typename OutputImageType::Pointer p_clOutputImage = p_clDistanceMapFilter->GetOutput();
  typename InputImageType::Pointer p_clVoronoiImage = p_clDistanceMapFilter->GetVoronoiMap();

  PairType clPair;
  clPair.first = p_clOutputImage;
  clPair.second = p_clVoronoiImage;

  return clPair;
}

template<typename PixelType>
bool GetSlice(typename itk::Image<PixelType, 2>::Pointer p_clSlice, typename itk::Image<PixelType, 3>::Pointer p_clVolume, unsigned int z) {
  typedef itk::Image<PixelType, 2> ImageType2D;
  typedef itk::Image<PixelType, 3> ImageType3D;
  typedef typename ImageType2D::SizeType SizeType2D;
  typedef typename ImageType3D::SizeType SizeType3D;
  typedef typename ImageType2D::SpacingType SpacingType2D;
  typedef typename ImageType3D::SpacingType SpacingType3D;
  typedef typename ImageType2D::PointType PointType2D;
  typedef typename ImageType3D::PointType PointType3D;

  if (!p_clVolume || !p_clSlice)
    return false;

  if (z >= p_clVolume->GetBufferedRegion().GetSize()[2])
    return false;

  SizeType2D clSize2D;
  clSize2D[0] = p_clVolume->GetBufferedRegion().GetSize()[0];
  clSize2D[1] = p_clVolume->GetBufferedRegion().GetSize()[1];

  SpacingType2D clSpacing2D;
  clSpacing2D[0] = p_clVolume->GetSpacing()[0];
  clSpacing2D[1] = p_clVolume->GetSpacing()[1];

  PointType2D clOrigin2D;
  clOrigin2D[0] = p_clVolume->GetOrigin()[0];
  clOrigin2D[1] = p_clVolume->GetOrigin()[1];

  p_clSlice->SetRegions(clSize2D);
  p_clSlice->SetSpacing(clSpacing2D);
  p_clSlice->SetOrigin(clOrigin2D);

  PixelType * const p_clVolumePointer = p_clVolume->GetPixelContainer()->GetBufferPointer();

  if (p_clVolumePointer == NULL)
    return false;

  p_clSlice->GetPixelContainer()->SetImportPointer(p_clVolumePointer + z * clSize2D[0]*clSize2D[1], clSize2D[0]*clSize2D[1], false);

  return true;
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer LoadVolume(const char *p_cFileName) {
  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer p_clReader = ReaderType::New();

  if (!p_clReader) {
    std::cerr << "Error: Could not allocate ReaderType." << std::endl;
    return typename ImageType::Pointer();
  }

  p_clReader->SetFileName(p_cFileName);

  try {
    p_clReader->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  typename ImageType::Pointer p_clOutputImage = p_clReader->GetOutput();

  return p_clOutputImage;
}

template<typename PixelType>
bool SaveVolume(typename itk::Image<PixelType, 3>::Pointer p_clVolume, const char *p_cFileName, bool bCompress) {
  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  if (!p_clVolume || p_cFileName == NULL)
    return false;

  typename WriterType::Pointer p_clWriter = WriterType::New();

  if (!p_clWriter) {  
    std::cerr << "Error: Could not allocate WriterType." << std::endl;
    return false;
  }

  p_clWriter->SetFileName(p_cFileName);
  p_clWriter->SetInput(p_clVolume.GetPointer());

  p_clWriter->SetUseCompression(bCompress);

  try {
    p_clWriter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return false;
  }

  return true;
}

template<typename PixelType>
bool SaveSlice(typename itk::Image<PixelType, 2>::Pointer p_clVolume, const char *p_cFileName, bool bCompress) {
  typedef itk::Image<PixelType, 2> ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  if (!p_clVolume || p_cFileName == NULL)
    return false;

  typename WriterType::Pointer p_clWriter = WriterType::New();

  if (!p_clWriter) {  
    std::cerr << "Error: Could not allocate WriterType." << std::endl;
    return false;
  }

  p_clWriter->SetFileName(p_cFileName);
  p_clWriter->SetInput(p_clVolume.GetPointer());
  p_clWriter->SetUseCompression(bCompress);

  try {
    p_clWriter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return false;
  }

  return true;
}

template<typename PixelType>
typename itk::Image<PixelType, 2>::Pointer ResampleSlice(typename itk::Image<PixelType, 2>::Pointer p_clSlice, int iNewWidth, int iNewHeight, bool bNearest) {
  typedef itk::Image<PixelType, 2> ImageType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef itk::LinearInterpolateImageFunction<ImageType> LinearInterpolationType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NearestNeighborInterpolationType;

  if (!p_clSlice)
    return typename ImageType::Pointer();

  const SizeType &clOriginalSize = p_clSlice->GetBufferedRegion().GetSize();
  const SpacingType &clOriginalSpacing = p_clSlice->GetSpacing();

  if (clOriginalSize[0] == iNewWidth && clOriginalSize[1] == iNewHeight)
    return p_clSlice; // Nothing to do

  SizeType clNewSize = { iNewWidth, iNewHeight };
  SpacingType clNewSpacing;

  clNewSpacing[0] = (clOriginalSpacing[0] * clOriginalSize[0])/clNewSize[0];
  clNewSpacing[1] = (clOriginalSpacing[1] * clOriginalSize[1])/clNewSize[1];
  
  typedef itk::IdentityTransform<double, 2> TransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;

  typename ResampleImageFilterType::Pointer p_clResampler = ResampleImageFilterType::New();

  p_clResampler->SetInput(p_clSlice);
  p_clResampler->SetSize(clNewSize);
  p_clResampler->SetOutputSpacing(clNewSpacing);
  p_clResampler->SetOutputOrigin(p_clSlice->GetOrigin());
  p_clResampler->SetTransform(TransformType::New());

  if (bNearest)
    p_clResampler->SetInterpolator(NearestNeighborInterpolationType::New());
  else
    p_clResampler->SetInterpolator(LinearInterpolationType::New());

  p_clResampler->UpdateLargestPossibleRegion();

  typename ImageType::Pointer p_clNewSlice = p_clResampler->GetOutput();

  p_clNewSlice->SetOrigin(p_clSlice->GetOrigin());
  p_clNewSlice->SetMetaDataDictionary(p_clSlice->GetMetaDataDictionary());

  return p_clNewSlice;
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer ResampleSlices(typename itk::Image<PixelType, 3>::Pointer p_clVolume, int iNewWidth, int iNewHeight, bool bNearest) {
  typedef itk::Image<PixelType, 3> ImageType3D;
  typedef itk::Image<PixelType, 2> ImageType2D;
  typedef typename ImageType3D::SizeType SizeType3D;
  typedef typename ImageType3D::SpacingType SpacingType3D;

  if (!p_clVolume)
    return typename ImageType3D::Pointer();

  const SizeType3D &clOriginalSize = p_clVolume->GetBufferedRegion().GetSize();
  const SpacingType3D &clOriginalSpacing = p_clVolume->GetSpacing();

  if (clOriginalSize[0] == iNewWidth && clOriginalSize[1] == iNewHeight)
    return p_clVolume; // Nothing to do

  SizeType3D clNewSize = { iNewWidth, iNewHeight, clOriginalSize[2] };
  SpacingType3D clNewSpacing;

  clNewSpacing[0] = (clOriginalSpacing[0] * clOriginalSize[0])/clNewSize[0];
  clNewSpacing[1] = (clOriginalSpacing[1] * clOriginalSize[1])/clNewSize[1];
  clNewSpacing[2] = clOriginalSpacing[2];

  typename ImageType3D::Pointer p_clNewVolume = ImageType3D::New();


  p_clNewVolume->SetRegions(clNewSize);
  p_clNewVolume->SetSpacing(clNewSpacing);
  p_clNewVolume->SetOrigin(p_clVolume->GetOrigin());
  p_clNewVolume->SetMetaDataDictionary(p_clVolume->GetMetaDataDictionary());

  p_clNewVolume->Allocate();

  typename ImageType2D::Pointer p_clSlice = ImageType2D::New();

  const size_t outputSliceSize = iNewWidth * iNewHeight;

  for (int z = 0; z < clOriginalSize[2]; ++z) {
    GetSlice<PixelType>(p_clSlice, p_clVolume, z);

    typename ImageType2D::Pointer p_clNewSlice = ResampleSlice<PixelType>(p_clSlice, iNewWidth, iNewHeight, bNearest);

    PixelType * const p_destBuffer = p_clNewVolume->GetPixelContainer()->GetBufferPointer() + outputSliceSize * z;
    const PixelType * const p_srcBuffer = p_clNewSlice->GetPixelContainer()->GetBufferPointer();

    std::copy(p_srcBuffer, p_srcBuffer + outputSliceSize, p_destBuffer);
  }

  return p_clNewVolume;
}

template<typename PixelType>
typename itk::Image<PixelType, 3>::Pointer ResampleVolume(typename itk::Image<PixelType, 3>::Pointer p_clVolume, int iNewWidth, int iNewHeight, int iNewDepth, bool bNearest) {
  typedef itk::Image<PixelType, 3> ImageType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef itk::LinearInterpolateImageFunction<ImageType> LinearInterpolationType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NearestNeighborInterpolationType;

  if (!p_clVolume)
    return typename ImageType::Pointer();

  const SizeType &clOriginalSize = p_clVolume->GetBufferedRegion().GetSize();
  const SpacingType &clOriginalSpacing = p_clVolume->GetSpacing();

  if (clOriginalSize[0] == iNewWidth && clOriginalSize[1] == iNewHeight && clOriginalSize[2] == iNewDepth)
    return p_clVolume; // Nothing to do

  SizeType clNewSize = { iNewWidth, iNewHeight, iNewDepth };
  SpacingType clNewSpacing;

  clNewSpacing[0] = (clOriginalSpacing[0] * clOriginalSize[0])/clNewSize[0];
  clNewSpacing[1] = (clOriginalSpacing[1] * clOriginalSize[1])/clNewSize[1];
  clNewSpacing[2] = (clOriginalSpacing[2] * clOriginalSize[2])/clNewSize[2];
  
  typedef itk::IdentityTransform<double, 3> TransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;

  typename ResampleImageFilterType::Pointer p_clResampler = ResampleImageFilterType::New();

  p_clResampler->SetInput(p_clVolume);
  p_clResampler->SetSize(clNewSize);
  p_clResampler->SetOutputSpacing(clNewSpacing);
  p_clResampler->SetOutputOrigin(p_clVolume->GetOrigin());
  p_clResampler->SetTransform(TransformType::New());

  if (bNearest)
    p_clResampler->SetInterpolator(NearestNeighborInterpolationType::New());
  else
    p_clResampler->SetInterpolator(LinearInterpolationType::New());

  p_clResampler->UpdateLargestPossibleRegion();

  typename ImageType::Pointer p_clNewVolume = p_clResampler->GetOutput();

  p_clNewVolume->SetOrigin(p_clVolume->GetOrigin());
  p_clNewVolume->SetMetaDataDictionary(p_clVolume->GetMetaDataDictionary());

  return p_clNewVolume;
}

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadDicomImage(const char *p_cPath, const char *p_cSeries) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef typename ImageType::PixelContainer PixelContainerType;
  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;

  if (p_cPath == NULL) {
    std::cerr << "Error: NULL path passed." << std::endl;
    return typename ImageType::Pointer();
  }

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();
  if (!p_clImageIO) {
    std::cerr << "Error: Could not allocate ImageIOType." << std::endl;
    return typename ImageType::Pointer();
  }

  // 2D slice
  if (Dimension <= 2) {
    if (p_cSeries != NULL)
      std::cout << "Warning: Ignoring series information when loading DICOM slice." << std::endl;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    
    typename ReaderType::Pointer p_clReader = ReaderType::New();

    if (!p_clReader) {
      std::cerr << "Error: Could not allocate ReaderType." << std::endl;
      return typename ImageType::Pointer();
    }

    p_clReader->SetImageIO(p_clImageIO);
    p_clReader->SetFileName(p_cPath);

    try {
      p_clReader->Update();
    }
    catch (itk::ExceptionObject &e) {
      std::cerr << "Error: " << e << std::endl;
      return typename ImageType::Pointer();
    }

    typename ImageType::Pointer p_clOutputImage = p_clReader->GetOutput();
    p_clOutputImage->SetMetaDataDictionary(p_clImageIO->GetMetaDataDictionary());

    return p_clOutputImage;
  }

  // >= 3D volume
  if (!IsFolder(p_cPath)) {
    std::cout << "Info: Got file. Trying to load '" << p_cPath << "' as a DICOM slice." << std::endl;

    typedef itk::Image<PixelType, 2> SliceType;
    typename SliceType::Pointer p_clSlice = LoadDicomImage<PixelType, 2>(p_cPath, p_cSeries);

    if (!p_clSlice) {
      std::cerr << "Error: Could not load file." << std::endl;
      return typename ImageType::Pointer();
    }

    const itk::MetaDataDictionary &clDicomTags = p_clSlice->GetMetaDataDictionary();

    std::string strSeriesUID;
    if (!itk::ExposeMetaData<std::string>(clDicomTags, "0020|000e", strSeriesUID)) {
      std::cerr << "Error: Could not get series UID." << std::endl;
      return typename ImageType::Pointer();
    }

    std::cout << "Info: Found series UID: " << strSeriesUID << std::endl;

    return LoadDicomImage<PixelType, Dimension>(DirName(p_cPath).c_str(), strSeriesUID.c_str());
  }

  typename ReaderType::Pointer p_clReader = ReaderType::New();
  if (!p_clReader) {
    std::cerr << "Error: Could not allocate ReaderType." << std::endl;
    return typename ImageType::Pointer();
  }

  p_clReader->SetImageIO(p_clImageIO);

  NamesGeneratorType::Pointer p_clNamesGenerator = NamesGeneratorType::New();

  if (!p_clNamesGenerator) {
    std::cerr << "Error: Could not allocate NamesGeneratorType." << std::endl;
    return typename ImageType::Pointer();
  }

  // Use all the DICOM information to determine if a series is part of the volume (must be called before SetDirectory())
  p_clNamesGenerator->SetUseSeriesDetails(false);
  p_clNamesGenerator->SetDirectory(p_cPath);

  std::string strSeriesUID;

  if (p_cSeries == NULL) {
    try {
      const std::vector<std::string> &vSeries = p_clNamesGenerator->GetSeriesUIDs();

      if (vSeries.empty()) {
        std::cerr << "Error: No series found in the DICOM files." << std::endl;
        return typename ImageType::Pointer();
      }

      //std::cout << "Info: There are " << vSeries.size() << " series." << std::endl;

      strSeriesUID = vSeries.front();
    }
    catch (itk::ExceptionObject &e) {
      std::cerr << "Error: " << e << std::endl;
      return typename ImageType::Pointer();
    }
  }
  else {
    strSeriesUID = p_cSeries;
  }

  //std::cout << "Info: Extracting series " << strSeriesUID << " ..." << std::endl;

  try {
    std::vector<std::string> vFileNames = p_clNamesGenerator->GetFileNames(strSeriesUID);

    //std::cout << "Info: There are " << vFileNames.size() << " files." << std::endl;

    p_clReader->SetFileNames(vFileNames);
    p_clReader->Update();

    typename ImageType::Pointer p_clOutputImage = p_clReader->GetOutput();
    p_clOutputImage->SetMetaDataDictionary(p_clImageIO->GetMetaDataDictionary());

    return p_clOutputImage;
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  return typename ImageType::Pointer(); // Not reached
}

template<typename PixelType, unsigned int Dimension>
typename itk::Image<PixelType, Dimension>::Pointer LoadDicomImageByIndex(const char *p_cPattern, int iBeginIndex, int iEndIndex) {
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef typename ImageType::PixelContainer PixelContainerType;
  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::NumericSeriesFileNames NamesGeneratorType;

  if (strstr(p_cPattern, "%d") == NULL && iEndIndex < iBeginIndex)
    iEndIndex = iBeginIndex;

  if (iEndIndex < iBeginIndex) {
    char a_cFileNameBuffer[2048];

    iEndIndex = iBeginIndex;
    do {
      ++iEndIndex;

#ifdef _WIN32
      sprintf_s(a_cFileNameBuffer, p_cPattern, iEndIndex);
#else // !_WIN32
      snprintf(a_cFileNameBuffer, sizeof(a_cFileNameBuffer), p_cPattern, iEndIndex);
#endif // _WIN32
    } while(FileExists(a_cFileNameBuffer));

    --iEndIndex;
  }

  //std::cout << "Info: Extracting series " << p_cPattern << ": begin = " << iBeginIndex << ", end = " << iEndIndex << " ..." << std::endl;

  NamesGeneratorType::Pointer p_clNamesGenerator = NamesGeneratorType::New();

  if (!p_clNamesGenerator) {
    std::cerr << "Error: Could not allocate NamesGeneratorType." << std::endl;
    return typename ImageType::Pointer();
  }

  p_clNamesGenerator->SetStartIndex(iBeginIndex);
  p_clNamesGenerator->SetEndIndex(iEndIndex);
  p_clNamesGenerator->SetSeriesFormat(p_cPattern);

  std::vector<std::string> vFileNames = p_clNamesGenerator->GetFileNames();
  //std::cout << "Info: There are " << vFileNames.size() << " files." << std::endl;

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();

  if (!p_clImageIO) {
    std::cerr << "Error: Could not allocate ImageIOType." << std::endl;
    return typename ImageType::Pointer();
  }

  typename ReaderType::Pointer p_clReader = ReaderType::New();

  if (!p_clReader) {
    std::cerr << "Error: Could not allocate ReaderType." << std::endl;
    return typename ImageType::Pointer();
  }

  p_clReader->SetImageIO(p_clImageIO);
  p_clReader->SetFileNames(vFileNames);

  try {
    p_clReader->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: Could not read series." << std::endl;
    std::cerr << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  return p_clReader->GetOutput();
}

// Work around broken 4D support in ITK
std::vector<std::string> GetTemporalFileNames(const char *p_cPath, const char *p_cSeries);

template<typename PixelType>
typename itk::Image<PixelType, 4>::Pointer LoadDicomVideo(const char *p_cPath, const char *p_cSeries) {
  typedef itk::Image<PixelType, 4> ImageType;
  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;

  if (p_cPath == NULL) {
    std::cerr << "Error: NULL path passed." << std::endl;
    return typename ImageType::Pointer();
  }

  if (!IsFolder(p_cPath)) {
    std::cout << "Info: Got file. Trying to load '" << p_cPath << "' as a DICOM slice." << std::endl;

    typedef itk::Image<PixelType, 2> SliceType;
    typename SliceType::Pointer p_clSlice = LoadDicomImage<PixelType, 2>(p_cPath, p_cSeries);

    if (!p_clSlice) {
      std::cerr << "Error: Could not load file." << std::endl;
      return typename ImageType::Pointer();
    }

    const itk::MetaDataDictionary &clDicomTags = p_clSlice->GetMetaDataDictionary();

    std::string strSeriesUID;
    if (!itk::ExposeMetaData<std::string>(clDicomTags, "0020|000e", strSeriesUID)) {
      std::cerr << "Error: Could not get series UID." << std::endl;
      return typename ImageType::Pointer();
    }

    std::cout << "Info: Found series UID: " << strSeriesUID << std::endl;

    return LoadDicomVideo<PixelType>(DirName(p_cPath).c_str(), strSeriesUID.c_str()); 
  }

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();
  if (!p_clImageIO) {
    std::cerr << "Error: Could not allocate ImageIOType." << std::endl;
    return typename ImageType::Pointer();
  }

  typename ReaderType::Pointer p_clReader = ReaderType::New();
  if (!p_clReader) {
    std::cerr << "Error: Could not allocate ReaderType." << std::endl;
    return typename ImageType::Pointer();
  }

  p_clReader->SetImageIO(p_clImageIO);

  std::vector<std::string> vFileNames = GetTemporalFileNames(p_cPath, p_cSeries);

  //for (size_t i = 0; i < vFileNames.size(); ++i) {
  //  std::cout << i << ": '" << vFileNames[i] << "'" << std::endl;
  //}

  typename ImageType::Pointer p_clVideo;

  try {
    //std::cout << "Info: There are " << vFileNames.size() << " files." << std::endl;

    p_clReader->SetFileNames(vFileNames);
    p_clReader->Update();

    p_clVideo = p_clReader->GetOutput();
    p_clVideo->SetMetaDataDictionary(p_clImageIO->GetMetaDataDictionary());
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return typename ImageType::Pointer();
  }

  typename ImageType::SizeType clSize = p_clVideo->GetBufferedRegion().GetSize();
  typename ImageType::SpacingType clSpacing = p_clVideo->GetSpacing();

  const itk::MetaDataDictionary &clDicomTags = p_clVideo->GetMetaDataDictionary();

  std::string strNumberOfTemporalPositions;
  if (!itk::ExposeMetaData(clDicomTags, "0020|0105", strNumberOfTemporalPositions)) {
    std::cerr << "Error: Could not read number of temporal positions (0020|0105)." << std::endl;
    return typename ImageType::Pointer();
  }

  Trim(strNumberOfTemporalPositions);

  char *p = NULL;
  const int iNumberOfTemporalPositions = strtol(strNumberOfTemporalPositions.c_str(), &p, 10);

  if (strNumberOfTemporalPositions.empty() || (p != NULL && *p != '\0')) {
    std::cerr << "Error: Could not parse '" << strNumberOfTemporalPositions << "' as an integer." << std::endl;
    return typename ImageType::Pointer();
  }

  if (iNumberOfTemporalPositions <= 0) {
    std::cerr << "Error: Invalid number of temporal positions (not a video?)." << std::endl;
    return typename ImageType::Pointer();
  }

  if ((clSize[2] % iNumberOfTemporalPositions) != 0) {
    std::cerr << "Error: " << iNumberOfTemporalPositions << " does not divide the number of slices " << clSize[2] << std::endl;
    return typename ImageType::Pointer();
  }

  const unsigned int uiNumVolumesInVideo = iNumberOfTemporalPositions;
  const unsigned int uiZDimension = (unsigned int)(clSize[2] / uiNumVolumesInVideo);

  clSize[2] = uiZDimension;
  clSize[3] = uiNumVolumesInVideo;

  p_clVideo->SetRegions(clSize);

  std::string strTemporalResolution;
  clSpacing[3] = 1.0;

  if (!itk::ExposeMetaData(clDicomTags, "0020|0110", strTemporalResolution)) {
    std::cerr << "Warning: Could not get temporal resolution." << std::endl;
  }
  else {
    Trim(strTemporalResolution);

    p = NULL;
    const double dTmp = strtod(strTemporalResolution.c_str(), &p);

    if (strTemporalResolution.empty() || (p != NULL && *p != '\0')) {
      std::cerr << "Warning: Could not parse '" << strTemporalResolution << "' as a double." << std::endl;
    }
    else {
      clSpacing[3] = dTmp;
    }
  }

  p_clVideo->SetSpacing(clSpacing);

  return p_clVideo;
}

template<typename PixelType>
bool SaveDicomSliceT(typename itk::Image<PixelType, 2>::Pointer p_clImage, const char *p_cFileName, bool bCompress) {
  typedef itk::Image<PixelType, 2> ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef itk::GDCMImageIO ImageIOType;

  if (!p_clImage || p_cFileName == NULL)
    return false;

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();

  if (!p_clImageIO) {
    std::cerr << "Error: Could not allocate ImageIOType." << std::endl;
    return false;
  }

  typename WriterType::Pointer p_clWriter = WriterType::New();

  if (!p_clWriter) {  
    std::cerr << "Error: Could not allocate WriterType." << std::endl;
    return false;
  }

  //p_clImage->SetMetaDataDictionary(clDicomTags);

  p_clWriter->SetImageIO(p_clImageIO);
  p_clWriter->SetFileName(p_cFileName);
  p_clWriter->SetInput(p_clImage.GetPointer());

  p_clImageIO->KeepOriginalUIDOn();

  p_clWriter->SetUseCompression(bCompress);

  p_clImageIO->SetMetaDataDictionary(p_clImage->GetMetaDataDictionary());

  try {
    p_clWriter->Update();
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return false;
  }

  return true;
}

template<typename Tp>
bool CsvRead(const char *p_cFileName, std::vector<std::vector<Tp> > &vMatrix, char cDelim = ',') {
  vMatrix.clear();

  std::ifstream csvStream(p_cFileName);

  if (!csvStream)
    return false;

  std::stringstream lineStream, tokenStream;
  std::string strLine, strToken;
  std::vector<Tp> vRow;

  while (std::getline(csvStream, strLine)) {
    vRow.clear();

    size_t p = strLine.find('\r');

    if (p != std::string::npos)
      strLine.erase(p);

    lineStream.clear();
    lineStream.str(strLine);

    while (std::getline(lineStream, strToken, cDelim)) {
      tokenStream.clear();
      tokenStream.str(strToken);

      Tp value = Tp();

      if (!(tokenStream >> value))
        return false;

      vRow.push_back(value);
    }

    if (vMatrix.size() > 0 && vMatrix[0].size() != vRow.size())
      return false;

    vMatrix.push_back(vRow);
  }

  return vMatrix.size() > 0 && vMatrix[0].size() > 0;
}

template<unsigned int Dimension>
bool NormalizeB2000(typename itk::Image<short,Dimension>::Pointer p_clB2000Image, typename itk::Image<short,Dimension>::Pointer p_clB0Image) {
  typedef itk::Image<short, Dimension> ImageType;

  if (!p_clB2000Image || !p_clB0Image)
    return false;

  const double dScale = 1000.0;

  const typename ImageType::SizeType &clB2000Size = p_clB2000Image->GetBufferedRegion().GetSize();
  const typename ImageType::SizeType &clB0Size = p_clB0Image->GetBufferedRegion().GetSize();

  if (clB2000Size != clB0Size) {
    std::cerr << "Error: Dimension mismatch between B0 and B2000 images." << std::endl;
    return false;
  }

  short * const p_sB2000Buffer = p_clB2000Image->GetPixelContainer()->GetBufferPointer();
  const size_t length = p_clB2000Image->GetPixelContainer()->Size();
  const short * const p_sB0Buffer = p_clB0Image->GetPixelContainer()->GetBufferPointer();

  for (size_t i = 0; i < length; ++i) {
    double dNewValue = 0.0;

    if (p_sB0Buffer[i] != 0)
      dNewValue = dScale * (double)p_sB2000Buffer[i] / (double)p_sB0Buffer[i] + 0.5;

    dNewValue = std::max(dNewValue, (double)std::numeric_limits<short>::min());
    dNewValue = std::min(dNewValue, (double)std::numeric_limits<short>::max());

    p_sB2000Buffer[i] = (short)dNewValue;
  }

  return true;
}

} // end namespace nih

#endif // COMMON_H
