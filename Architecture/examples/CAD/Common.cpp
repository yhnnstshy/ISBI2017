#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <direct.h>

#undef MoveFile
#undef CopyFile

#elif defined(__unix__)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <fts.h>
#include <glob.h>

#else // !_WIN32
#error "Unsupported"
#endif // _WIN32

#include "Common.h"

namespace nih {

namespace {
  std::string g_strLogFile = "console";
} // end anonymous namespace

int64_t Rand64() {
  // NOTE: Should use C++11 random features in the future
  // XXX: This may not actually be uniform!!!

#ifdef _WIN32
  // Windows is dumb. It only supports 15 bit random numbers from rand()
  return (int64_t)rand() |
        ((int64_t)rand() << 15) |
        ((int64_t)rand() << 30) |
        ((int64_t)rand() << 45);
#else // !_WIN32
  // Every other OS at least returns 31 bit random numbers
  return (int64_t)rand() |
        ((int64_t)rand() << 31);
#endif // _WIN32
}

bool MakeFolder(const char *p_cPath) {
#ifdef _WIN32
  return _mkdir(p_cPath) == 0;
#else // !_WIN32
  return mkdir(p_cPath, 0777) == 0;
#endif // _WIN32
}

bool DeleteFolder(const char *p_cPath) {
#ifdef _WIN32
  return _rmdir(p_cPath) == 0;
#else // !_WIN32
  return rmdir(p_cPath) == 0;
#endif // _WIN32
}

#ifdef _WIN32
bool FileExists(const char *p_cPath) {
  return GetFileAttributes(p_cPath) != INVALID_FILE_ATTRIBUTES;
}
#else // !_WIN32
bool FileExists(const char *p_cPath) {
  struct stat buf;
  memset(&buf, 0, sizeof(buf));
  return stat(p_cPath, &buf) == 0;
}
#endif

bool IsFolder(const char *p_cPath) {
#ifdef _WIN32
  const DWORD dwFileAttr = GetFileAttributes(p_cPath);
  return (dwFileAttr != INVALID_FILE_ATTRIBUTES) && (dwFileAttr & FILE_ATTRIBUTE_DIRECTORY);
#else // !_WIN32
  struct stat buf;
  memset(&buf, 0, sizeof(buf));
  return stat(p_cPath, &buf) == 0 && S_ISDIR(buf.st_mode);  
#endif // _WIN32
}

bool MoveFile(const char *p_cSource, const char *p_cDest) {
#ifdef _WIN32
  return ::MoveFileExA(p_cSource, p_cDest, MOVEFILE_COPY_ALLOWED | MOVEFILE_WRITE_THROUGH | MOVEFILE_FAIL_IF_NOT_TRACKABLE) != 0;
#else // !_WIN32
  return rename(p_cSource, p_cDest) == 0;
#endif // _WIN32
}

#ifdef _WIN32
bool CopyFile(const char *p_cSource, const char *p_cDest) {
  return ::CopyFileA(p_cSource, p_cDest, FALSE);
}
#else // !_WIN32

bool CopyFile(const char *p_cSource, const char *p_cDest) {
  enum { COPY_BUFFER_SIZE = 1024*1024 };

  int iFromFd = open(p_cSource, O_RDONLY, 0);

  if (iFromFd == -1)
    return false;

  struct stat stBuf;
  memset(&stBuf, 0, sizeof(stBuf));

  if (fstat(iFromFd, &stBuf) != 0) {
    close(iFromFd);
    return false;
  }

  int iToFd = open(p_cDest, O_WRONLY | O_TRUNC | O_CREAT, stBuf.st_mode & ~(S_ISUID | S_ISGID));

  if (iToFd == -1) {
    close(iFromFd);
    return false;
  }

  bool bSucceed = true;
  char a_cBuffer[COPY_BUFFER_SIZE];
  ssize_t readSize = 0;
  while ((readSize = read(iFromFd, a_cBuffer, sizeof(a_cBuffer))) > 0) {
    ssize_t writeSize = 0;
    char *p_cPos = a_cBuffer;

    while ((writeSize = write(iToFd, p_cPos, readSize)) > 0) {
      if (writeSize < 0) {
        bSucceed = false;
        break;
      }

      readSize -= writeSize;
      p_cPos += writeSize;
    }

    if (!bSucceed)
      break;
  }

  close(iToFd);
  close(iFromFd);

  return bSucceed;
}

#endif // _WIN32

#ifdef _WIN32
void FindDicomFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder) {
  WIN32_FIND_DATA clData;

  std::string strPattern = p_cRootFolder;
  strPattern += "\\*";

  HANDLE hFindData = FindFirstFile(strPattern.c_str(), &clData);

  if (hFindData == NULL) {
    std::cerr << "Error: Could not walk '" << p_cRootFolder << "'." << std::endl;
    return;
  }

  bool bFoundSlice = false;

  do {
    if (strcmp(clData.cFileName, ".") == 0 || strcmp(clData.cFileName, "..") == 0)
      continue;

    const bool bIsFolder = (clData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
    std::string strFullPath = p_cRootFolder;
    strFullPath += '\\';
    strFullPath += clData.cFileName;

    //std::cout << "Found: " << strFullPath << std::endl;

    if (bIsFolder) {
      // Recurse
      FindDicomFolders(vFolders, strFullPath.c_str());
      continue;
    }

    // Already found a slice for this folder?
    if (bFoundSlice)
      continue;

    // Otherwise it's a file
    if (LoadDicomSlice(strFullPath.c_str())) {
      vFolders.push_back(p_cRootFolder);
      bFoundSlice = true;
    }

  } while (FindNextFile(hFindData, &clData));

  FindClose(hFindData);
}
#else // !_WIN32

void FindDicomFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder) {
  const char * const p_cPathArgv[2] = { p_cRootFolder, NULL };
  FTS *p_ftsCtx = fts_open((char * const *)p_cPathArgv, FTS_COMFOLLOW | FTS_NOCHDIR | FTS_LOGICAL | FTS_NOSTAT, NULL);

  if (p_ftsCtx == NULL)
    return;

  FTSENT *p_ftsEnt = NULL;
  while ((p_ftsEnt = fts_read(p_ftsCtx)) != NULL) {
    switch (p_ftsEnt->fts_info) {
    case FTS_F:
      if (p_ftsEnt->fts_parent->fts_number != 0)
        break;

      if (LoadDicomSlice(p_ftsEnt->fts_path)) {
        vFolders.push_back(p_ftsEnt->fts_parent->fts_path);
        p_ftsEnt->fts_parent->fts_number = 1;
      }

      break;
    }
  }

  fts_close(p_ftsCtx);
}

#endif // _WIN32

#ifdef _WIN32

void FindFiles(std::vector<std::string> &vFiles, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive) {
  WIN32_FIND_DATA clData;

  std::string strPattern = p_cRootFolder;
  strPattern += '\\';
  strPattern += p_cPattern;

  HANDLE hFindData = FindFirstFile(strPattern.c_str(), &clData);

  if (hFindData == NULL) {
    std::cerr << "Error: Could not walk '" << p_cRootFolder << "'." << std::endl;
    return;
  }

  do {
    if (strcmp(clData.cFileName, ".") == 0 || strcmp(clData.cFileName, "..") == 0)
      continue;

    const bool bIsFolder = (clData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
    std::string strFullPath = p_cRootFolder;
    strFullPath += '\\';
    strFullPath += clData.cFileName;

    //std::cout << "Found: " << strFullPath << std::endl;

    if (!bIsFolder)
      vFiles.push_back(strFullPath);

  } while (FindNextFile(hFindData, &clData));

  FindClose(hFindData);

  if (bRecursive) {
    strPattern = p_cRootFolder;
    strPattern += "\\*";

    hFindData = FindFirstFile(strPattern.c_str(), &clData);
    
    do {
      if (strcmp(clData.cFileName, ".") == 0 || strcmp(clData.cFileName, "..") == 0)
        continue;

      const bool bIsFolder = (clData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
      std::string strFullPath = p_cRootFolder;
      strFullPath += '\\';
      strFullPath += clData.cFileName;

      if (bIsFolder)
        FindFiles(vFiles, strFullPath.c_str(), p_cPattern, bRecursive);

    } while (FindNextFile(hFindData, &clData));

    FindClose(hFindData);
  }
}

#else

void FindFiles(std::vector<std::string> &vFiles, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive) {
  glob_t globCtx;

  std::string strPattern = p_cRootFolder;
  strPattern += '/';
  strPattern += p_cPattern;

  memset(&globCtx, 0, sizeof(globCtx));
  if (glob(strPattern.c_str(), GLOB_TILDE, NULL, &globCtx) == 0) {
    for (size_t i = 0; i < globCtx.gl_pathc; ++i) {
      if (!IsFolder(globCtx.gl_pathv[i]))
        vFiles.push_back(globCtx.gl_pathv[i]);
    }
  }

  globfree(&globCtx);

  if (bRecursive) {
    memset(&globCtx, 0, sizeof(globCtx));

    strPattern = p_cRootFolder;
    strPattern += "/*";

    if (glob(strPattern.c_str(), GLOB_TILDE, NULL, &globCtx) == 0) {
      for (size_t i = 0; i < globCtx.gl_pathc; ++i) {
        if (strcmp(globCtx.gl_pathv[i], ".") != 0 && strcmp(globCtx.gl_pathv[i], "..") != 0 && IsFolder(globCtx.gl_pathv[i]))
          FindFiles(vFiles, globCtx.gl_pathv[i], p_cPattern, bRecursive);
      }
    }

    globfree(&globCtx);
  }
}

#endif // _WIN32

#ifdef _WIN32
void FindFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive) {
  WIN32_FIND_DATA clData;

  std::string strPattern = p_cRootFolder;
  strPattern += '\\';
  strPattern += p_cPattern;

  HANDLE hFindData = FindFirstFile(strPattern.c_str(), &clData);

  if (hFindData == NULL) {
    std::cerr << "Error: Could not walk '" << p_cRootFolder << "'." << std::endl;
    return;
  }

  do {
    if (strcmp(clData.cFileName, ".") == 0 || strcmp(clData.cFileName, "..") == 0)
      continue;

    const bool bIsFolder = (clData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
    std::string strFullPath = p_cRootFolder;
    strFullPath += '\\';
    strFullPath += clData.cFileName;

    //std::cout << "Found: " << strFullPath << std::endl;

    if (bIsFolder)
      vFolders.push_back(strFullPath);

  } while (FindNextFile(hFindData, &clData));

  FindClose(hFindData);

  if (bRecursive) {
    strPattern = p_cRootFolder;
    strPattern += "\\*";

    hFindData = FindFirstFile(strPattern.c_str(), &clData);
    
    do {
      if (strcmp(clData.cFileName, ".") == 0 || strcmp(clData.cFileName, "..") == 0)
        continue;

      const bool bIsFolder = (clData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
      std::string strFullPath = p_cRootFolder;
      strFullPath += '\\';
      strFullPath += clData.cFileName;

      if (bIsFolder)
        FindFiles(vFolders, strFullPath.c_str(), p_cPattern, bRecursive);

    } while (FindNextFile(hFindData, &clData));

    FindClose(hFindData);
  }
}
#else // !_WIN32

void FindFolders(std::vector<std::string> &vFolders, const char *p_cRootFolder, const char *p_cPattern, bool bRecursive) {
  glob_t globCtx;

  std::string strPattern = p_cRootFolder;
  strPattern += '/';
  strPattern += p_cPattern;

  memset(&globCtx, 0, sizeof(globCtx));
  if (glob(strPattern.c_str(), GLOB_TILDE, NULL, &globCtx) == 0) {
    for (size_t i = 0; i < globCtx.gl_pathc; ++i) {
      if (IsFolder(globCtx.gl_pathv[i]))
        vFolders.push_back(globCtx.gl_pathv[i]);
    }
  }

  globfree(&globCtx);

  if (bRecursive) {
    memset(&globCtx, 0, sizeof(globCtx));

    strPattern = p_cRootFolder;
    strPattern += "/*";

    if (glob(strPattern.c_str(), GLOB_TILDE, NULL, &globCtx) == 0) {
      for (size_t i = 0; i < globCtx.gl_pathc; ++i) {
        if (strcmp(globCtx.gl_pathv[i], ".") != 0 && strcmp(globCtx.gl_pathv[i], "..") != 0 && IsFolder(globCtx.gl_pathv[i]))
          FindFolders(vFolders, globCtx.gl_pathv[i], p_cPattern, bRecursive);
      }
    }

    globfree(&globCtx);
  }
}

#endif // _WIN32

std::string GetCurrentTimeString() {
  char a_cTimeBuf[128] = "";
  time_t tTime = 0;

  time(&tTime);

// NOTE: We need to do this in a thread-safe portable way (ctime() is NOT thread safe)

#ifdef _WIN32
  ctime_s(a_cTimeBuf, sizeof(a_cTimeBuf), &tTime);
#else // !_WIN32
  // NOTE: FreeBSD man page says that the character buffer should be at least 26 characters
  ctime_r(&tTime, a_cTimeBuf);
#endif // _WIN32

  return std::string(a_cTimeBuf);
}

std::string DirName(const std::string &strFilePath) {
  if (strFilePath.empty())
    return ".";

  std::string strNewFilePath = strFilePath;

#ifdef _WIN32
  while (strNewFilePath.size() > 1 && (strNewFilePath.back() == '/' || strNewFilePath.back() == '\\'))
    strNewFilePath.pop_back();

  if (strNewFilePath == "/" || strNewFilePath == "\\")
    return strNewFilePath;

  size_t p = strNewFilePath.find_last_of("/\\");
#else // !_WIN32
  while (strNewFilePath.size() > 1 && strNewFilePath.back() == '/')
    strNewFilePath.pop_back();

  if (strNewFilePath == "/")
    return strNewFilePath;

  size_t p = strNewFilePath.find_last_of("/");
#endif // _WIN32

  if (p == std::string::npos)
    return ".";

  strNewFilePath.erase(p);

  return strNewFilePath;
}

std::string BaseName(const std::string &strFilePath) {
  if (strFilePath.empty())
    return "";

  std::string strNewFilePath = strFilePath;

#ifdef _WIN32
  while (strNewFilePath.size() > 1 && (strNewFilePath.back() == '/' || strNewFilePath.back() == '\\'))
    strNewFilePath.pop_back();

  if (strNewFilePath == "/" || strNewFilePath == "\\")
    return strNewFilePath;

  size_t p = strNewFilePath.find_last_of("/\\");
#else // !_WIN32
  while (strNewFilePath.size() > 1 && strNewFilePath.back() == '/')
    strNewFilePath.pop_back();

  if (strNewFilePath == "/")
    return strNewFilePath;

  size_t p = strNewFilePath.find_last_of("/");
#endif // _WIN32

  if (p == std::string::npos)
    return strNewFilePath;

  strNewFilePath.erase(0,p+1);

  return strNewFilePath;
}

std::string StripExtension(const std::string &strFilePath) {


  // It could happen that the last . occurs before a path delimeter
  size_t p = strFilePath.rfind('.');

  if (p == std::string::npos)
    return strFilePath;

#ifdef _WIN32
  size_t q = strFilePath.find_last_of("/\\");
#else // !_WIN32
  size_t q = strFilePath.rfind('/');
#endif // _WIN32

  if (q != std::string::npos && p < q)
    return strFilePath;

  std::string strNewFilePath;

  strNewFilePath.assign(strFilePath, 0, p);

  return strNewFilePath;
}

void Trim(std::string &strString) {
  size_t p = strString.find_first_not_of(" \t\r\n");
  if (p != std::string::npos)
    strString.erase(0, p);

  p = strString.find_last_not_of(" \t\r\n");
  if (p != std::string::npos && p+1 < strString.size())
    strString.erase(p+1);
}

void SetLogFile(const char *p_cLogFile) { 
  if (p_cLogFile == NULL)
    g_strLogFile.clear();
  else
    g_strLogFile = p_cLogFile;
}

void ClearLogFile() {
  if (g_strLogFile.empty() || g_strLogFile == "console")
    return;

  // Truncate
  std::ofstream truncStream(g_strLogFile.c_str(), std::ofstream::trunc);
}

void Log(const char *p_cFormat, ...) {
  if (g_strLogFile.empty())
    return;

  va_list ap;

  if (g_strLogFile == "console") {
    va_start(ap, p_cFormat);    
    vprintf(p_cFormat, ap);
    va_end(ap);
  }
  else {
    FILE *pFile = fopen(g_strLogFile.c_str(), "a");

    if (pFile == NULL)
      return;

    va_start(ap, p_cFormat);
    vfprintf(pFile, p_cFormat, ap);
    va_end(ap);

    fclose(pFile);
  }
}

void RegisterImageFormats() {
  itk::MetaImageIOFactory::RegisterOneFactory();
  itk::GDCMImageIOFactory::RegisterOneFactory();
  itk::NiftiImageIOFactory::RegisterOneFactory();
  itk::NrrdImageIOFactory::RegisterOneFactory();
}

// Work around broken 4D support in ITK
std::vector<std::string> GetTemporalFileNames(const char *p_cPath, const char *p_cSeries) {
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();

  if (!p_clImageIO) {
    std::cerr << "Error: Failed to allocate ImageIOType." << std::endl;
    return std::vector<std::string>();
  }

  NamesGeneratorType::Pointer p_clNamesGenerator = NamesGeneratorType::New();

  if (!p_clNamesGenerator) {
    std::cerr << "Error: Failed to allocate NamesGeneratorType." << std::endl;
    return std::vector<std::string>();
  }

  p_clNamesGenerator->SetUseSeriesDetails(false);
  p_clNamesGenerator->SetDirectory(p_cPath);

  std::string strSeriesUID;

  if (p_cSeries == NULL) {
    try {
      std::vector<std::string> vAllSeriesUIDs = p_clNamesGenerator->GetSeriesUIDs();
      if (vAllSeriesUIDs.empty()) {
        std::cerr << "Error: No series found in '" << p_cPath << "'." << std::endl;
        return std::vector<std::string>();
      }

      strSeriesUID = vAllSeriesUIDs.front();
    }
    catch (itk::ExceptionObject &e) {
      std::cerr << "Error: " << e << std::endl;
      return std::vector<std::string>();
    }
  }
  else {
    strSeriesUID = p_cSeries;
  }

  // First get all files names for ths series

  std::unordered_set<std::string> sAllFiles;
  try {
    std::vector<std::string> vAllFiles = p_clNamesGenerator->GetFileNames(strSeriesUID);
    sAllFiles.insert(vAllFiles.begin(), vAllFiles.end());
  }
  catch (itk::ExceptionObject &e) {
    std::cout << "Error: " << e << std::endl;
    return std::vector<std::string>();
  }

  if (sAllFiles.empty()) {
    std::cerr << "Error: No files found for series '" << strSeriesUID << "'." << std::endl;
    return std::vector<std::string>();
  }

  // Map temporal ID
  typedef std::map<int, std::vector<std::string> > MapType;

  MapType mVolumeMap;

  // XXX: Create new names generator since the old data is re-used otherwise.
  p_clNamesGenerator = NamesGeneratorType::New();
  if (!p_clNamesGenerator) {
    std::cerr << "Error: Could not allocate NamesGeneratorType." << std::endl;
    return std::vector<std::string>();
  }

  p_clNamesGenerator->SetUseSeriesDetails(true);
  p_clNamesGenerator->AddSeriesRestriction("0020|0100");
  p_clNamesGenerator->SetDirectory(p_cPath);

  // NOTE: These are NOT the same series UID since they are appended with additional information!
  std::vector<std::string> vAllSeriesUIDs = p_clNamesGenerator->GetSeriesUIDs();

  std::cout << "There are " << vAllSeriesUIDs.size() << " series." << std::endl;

  for (size_t i = 0; i < vAllSeriesUIDs.size(); ++i) {
    const std::string &strSeriesUID = vAllSeriesUIDs[i];

    try {
      std::vector<std::string> vFileNames = p_clNamesGenerator->GetFileNames(strSeriesUID);
      if (vFileNames.empty()) {
        // XXX: Shouldn't happen!
        std::cerr << "Error: Somehow there are no files for the automatically detected series UID: '" << strSeriesUID << "'." << std::endl;
        return std::vector<std::string>();
      }

      if (sAllFiles.find(vFileNames[0]) != sAllFiles.end()) {
        itk::Image<short, 2>::Pointer p_clSlice = LoadDicomSlice(vFileNames[0].c_str());
        if (!p_clSlice) {
          // XXX: Shouldn't happen!
          std::cerr << "Error: Could not load slice from automatically detected filename: '" << vFileNames[0] << "'." << std::endl;
          return std::vector<std::string>();
        }

        // Get temporal position id
        std::string strTemporalPosition;
        if (!itk::ExposeMetaData(p_clSlice->GetMetaDataDictionary(), "0020|0100", strTemporalPosition)) {
          std::cerr << "Error: Could not get temporal position id for DICOM slice '" << vFileNames[0] << "'." << std::endl;
          return std::vector<std::string>();
        }

        Trim(strTemporalPosition);

        const int iTemporalPosition = strtol(strTemporalPosition.c_str(), NULL, 10); // XXX: Just going to assume this is correct!
        
        mVolumeMap[iTemporalPosition].swap(vFileNames); // Should only occur once!
      }
    }
    catch (itk::ExceptionObject &e) {
      // XXX: Shouldn't happen!
      std::cerr << "Error: " << e << std::endl;
      return std::vector<std::string>();
    }
  }

  std::vector<std::string> vOrderedFiles;

  for (MapType::iterator itr = mVolumeMap.begin(); itr != mVolumeMap.end(); ++itr) {
    const std::vector<std::string> &vFileNames = itr->second;
    vOrderedFiles.insert(vOrderedFiles.end(), vFileNames.begin(), vFileNames.end());
  }

  return vOrderedFiles;
}

bool GetBValueFileNames(std::map<double, std::vector<std::string> > &mBValueMap, const char *p_cPath) {
  typedef itk::Image<short, 2> ImageType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;

  NamesGeneratorType::Pointer p_clNamesGenerator = NamesGeneratorType::New();

  p_clNamesGenerator->SetUseSeriesDetails(true);
  p_clNamesGenerator->AddSeriesRestriction("0018|9087"); // B-value tag
  p_clNamesGenerator->SetDirectory(p_cPath);

  std::vector<std::string> vSeries;
  try {
    vSeries = p_clNamesGenerator->GetSeriesUIDs();

    if (vSeries.empty()) {
      std::cerr << "Error: No series found in the DICOM files." << std::endl;
      return false;
    }

    std::cout << "Info: There are " << vSeries.size() << " series." << std::endl;
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < vSeries.size(); ++i) {
    const std::string &strSeries = vSeries[i];
    std::vector<std::string> vFileNames = p_clNamesGenerator->GetFileNames(vSeries[i]);

    if (vFileNames.empty())
      continue; // Shouldn't happen

    ImageType::Pointer p_clImage = LoadDicomSlice(vFileNames[0].c_str());
    if (!p_clImage) {
      std::cerr << "Error: Could not read slice from series '" << strSeries << "'." << std::endl;
      continue;
    }

    const itk::MetaDataDictionary &clDicomTags = p_clImage->GetMetaDataDictionary();

    std::string strBValue;

    if (!itk::ExposeMetaData<std::string>(clDicomTags, "0018|9087", strBValue))
      continue;

    double dBValue = -1.0;
    std::stringstream valueStream;
    valueStream.str(strBValue);

    if (!(valueStream >> dBValue))
      continue;

    std::cout << "Info: B value " << dBValue << " = " << strSeries << std::endl;

    mBValueMap[dBValue] = vFileNames;
  }

  return mBValueMap.size() > 0;
}

itk::Image<short, 3>::Pointer LoadB2000Volume(const char *p_cPath, bool bNormalize) {
  // This series has all sequences stacked together

  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;

  if (!IsFolder(p_cPath))
    return LoadB2000Volume(DirName(p_cPath).c_str(), bNormalize);

  typedef std::map<double, std::vector<std::string> > BValueMapType;

  BValueMapType mBValueMap;
  if (!GetBValueFileNames(mBValueMap, p_cPath)) {
    std::cerr << "Error: Could not get B value files." << std::endl;
    return ImageType::Pointer();
  }

  std::vector<std::string> vB0Files, vB2000Files;

  for (BValueMapType::iterator itr = mBValueMap.begin(); itr != mBValueMap.end(); ++itr) {
    const double dBValue = itr->first;
    const std::vector<std::string> &vFileNames = itr->second;

    if (std::abs(dBValue - 2000.0) < 10)
      vB2000Files = vFileNames;
    else if (std::abs(dBValue) < 1e-1)
      vB0Files = vFileNames;
  }

  if (vB2000Files.empty()) {
    std::cerr << "Error: No B2000 files." << std::endl; 
    return ImageType::Pointer();
  }

  ReaderType::Pointer p_clReader = ReaderType::New();
  if (!p_clReader) {
    std::cerr << "Error: Could not allocate ReaderType." << std::endl;
    return ImageType::Pointer();
  }

  ImageIOType::Pointer p_clImageIO = ImageIOType::New();
  if (!p_clImageIO) {
    std::cerr << "Error: Could not allocate ImageIOType." << std::endl;
    return ImageType::Pointer();
  }

  p_clReader->SetImageIO(p_clImageIO);
  
  ImageType::Pointer p_clB2000Image, p_clB0Image;

  p_clReader->SetFileNames(vB2000Files);

  try {
    p_clReader->Update();
    p_clB2000Image = p_clReader->GetOutput();
    p_clB2000Image->SetMetaDataDictionary(p_clImageIO->GetMetaDataDictionary());
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: Could not read B2000 volume: " << e << std::endl;
    return ImageType::Pointer();
  }

  if (!bNormalize)
    return p_clB2000Image;

  if (vB0Files.empty()) {
    std::cerr << "Error: No B0 files." << std::endl;
    return ImageType::Pointer();
  }

  p_clReader->SetFileNames(vB0Files);

  try {
    p_clReader->Update();
    p_clB0Image = p_clReader->GetOutput();
    p_clB0Image->SetMetaDataDictionary(p_clImageIO->GetMetaDataDictionary());
  }
  catch (itk::ExceptionObject &e) {
    std::cerr << "Error: Could not read B0 volume: " << e << std::endl;
    return ImageType::Pointer();
  }

  if (!NormalizeB2000<3>(p_clB2000Image, p_clB0Image)) {
    std::cerr << "Error: Failed to normalize B2000 image." << std::endl;
    return ImageType::Pointer();
  }

  return p_clB2000Image;
}

} // end namespace nih
