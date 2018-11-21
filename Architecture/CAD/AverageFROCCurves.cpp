#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <limits>
#include <string>
#include <sstream>
#include "PiecewiseLinearInterpolation.h"
#include "bsdgetopt.h"

struct OperatingPoint {
  float fThreshold;
  float fTruePositiveRate;
  float fFalsePositiveCount;

  float fThresholdMin, fThresholdMax, fThresholdStd;
  float fTruePositiveRateMin, fTruePositiveRateMax, fTruePositiveRateStd;

  OperatingPoint() {
    fThreshold = fThresholdMin = fThresholdMax = fThresholdStd = 0.0f;
    fTruePositiveRate = fTruePositiveRateMin = fTruePositiveRateMax = fTruePositiveRateStd = 0.0f;
    fFalsePositiveCount = 0.0f;
  }

  OperatingPoint(float fThreshold_, float fTruePositiveRate_, float fFalsePositiveRate_) {
    fThreshold = fThresholdMin = fThresholdMax = fThreshold_;
    fTruePositiveRate = fTruePositiveRateMin = fTruePositiveRateMax = fTruePositiveRate_;
    fFalsePositiveCount = fFalsePositiveRate_;

    fTruePositiveRateStd = fThresholdStd = 0.0f;
  }

  bool operator<(const OperatingPoint &stOther) const {
    return fFalsePositiveCount < stOther.fFalsePositiveCount;
  }
};

std::vector<OperatingPoint> LoadFROC(const char *p_cFileName);
void SaveFROC(const char *p_cFileName, std::vector<OperatingPoint> vCurve);
std::vector<OperatingPoint> ResampleCurve(const std::vector<OperatingPoint> &vCurve, const std::set<float> &sFPCounts);
std::vector<OperatingPoint> AverageFROCCurves(const std::vector<std::vector<OperatingPoint> > &vCurves);
OperatingPoint InterpolateOverTPR(const std::vector<OperatingPoint> &vCurve, float fTPR);
OperatingPoint InterpolateOverFPC(const std::vector<OperatingPoint> &vCurve, float fFPC);
OperatingPoint InterpolateOverProb(const std::vector<OperatingPoint> &vCurve, float fProb);

void Usage(const char *p_cArg0) {
  std::cerr << "Usage: " << p_cArg0 << " [-t interpolateTPR] [-f interpolateFPCount] [-p probThreshold] [-h] [-o stdout] FROCCurve.txt [FROCCurve2.txt ...]" << std::endl;
  exit(1);
}

int main(int argc, char **argv) {
  const char * const p_cArg0 = argv[0];

  std::string strOutputFile = "stdout";

  float fTruePositiveRate = -1.0f;
  float fFalsePositiveCount = -1.0f;
  float fProb = -1.0f;

  int c = 0;
  while ((c = getopt(argc, argv, "f:ho:p:t:")) != -1) {
    switch (c) {
    case 'f':
      {
        char *p = NULL; 
        fFalsePositiveCount = (float)strtod(optarg, &p);
        if (*p != '\0') 
          Usage(p_cArg0);
      }
      break;
    case 'h':
      Usage(p_cArg0);
      break;
    case 'o':
      strOutputFile = optarg;
      break;
    case 'p':
      {
        char *p = NULL; 
        fProb = (float)strtod(optarg, &p);
        if (*p != '\0') 
          Usage(p_cArg0);
      }
      break;
    case 't':
      {
        char *p = NULL; 
        fTruePositiveRate = (float)strtod(optarg, &p);
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

  if (argc <= 0)
    Usage(p_cArg0);

  std::vector<std::vector<OperatingPoint> > vCurves(argc);

  for (int i = 0; i < argc; ++i) {
    const char * const p_cFileName = argv[i];

    vCurves[i] = LoadFROC(p_cFileName);
    if (vCurves[i].empty()) {
      std::cerr << "Error: Could not load file '" << p_cFileName << "'." << std::endl;
      return -1;
    }
  }

  std::vector<OperatingPoint> vMeanCurve = AverageFROCCurves(vCurves);

  if (fTruePositiveRate >= 0.0f) {
    OperatingPoint clResult = InterpolateOverTPR(vMeanCurve, fTruePositiveRate);
    std::cout << "True Positive Rate: " << fTruePositiveRate << " --> False Positive Count: " << clResult.fFalsePositiveCount << ", Threshold: " << clResult.fThreshold << std::endl;
    return 0;
  }

  if (fFalsePositiveCount >= 0.0f) {
    OperatingPoint clResult = InterpolateOverFPC(vMeanCurve, fFalsePositiveCount);
    std::cout << "False Positive Count: " << fFalsePositiveCount << " --> True Positive Rate: " << clResult.fTruePositiveRate << ", Threshold: " << clResult.fThreshold << std::endl;
    return 0;
  }

  if (fProb >= 0.0f) {
    OperatingPoint clResult = InterpolateOverProb(vMeanCurve, fProb);
    std::cout << "Threshold: " << fProb << " --> True Positive Rate: " << clResult.fTruePositiveRate << ", False Positive Count: " << clResult.fFalsePositiveCount << std::endl;
    return 0;
  }

  SaveFROC(strOutputFile.c_str(), vMeanCurve);

  return 0;
}

std::vector<OperatingPoint> LoadFROC(const char *p_cFileName) {
  std::ifstream fileStream(p_cFileName);

  if (!fileStream)
    return std::vector<OperatingPoint>();

  std::vector<OperatingPoint> vCurve;

  std::string strLine;

  std::stringstream lineStream;
  while (std::getline(fileStream, strLine)) {
    size_t p = strLine.find('#');
    if (p != std::string::npos)
      strLine.erase(p);

    lineStream.clear();
    lineStream.str(strLine);

    float fThreshold = 0.0f, fTruePositiveRate = 0.0f, fFalsePositiveRate = 0.0f;
    if (!(lineStream >> fThreshold >> fTruePositiveRate >> fFalsePositiveRate))
      continue;
    
    vCurve.push_back(OperatingPoint(fThreshold, fTruePositiveRate, fFalsePositiveRate));
  }

  return vCurve;
}

void SaveFROC(const char *p_cFileName, std::vector<OperatingPoint> vCurve) {
  std::ofstream fileStream;
  if (strcmp(p_cFileName, "stdout") != 0) {
    fileStream.open(p_cFileName, std::ofstream::out | std::ofstream::trunc);

    if (!fileStream.good()) {
      std::cerr << "Error: Could not open '" << p_cFileName << "'." << std::endl;
      return;
    }
  }

  std::ostream &os = fileStream.is_open() ? fileStream : std::cout;

  std::sort(vCurve.begin(), vCurve.end());

  os << "# Threshold\tTrue Positive Rate\tFalse Positive Count\tThresholdMin\tThresholdMax\tThresholdStd\tTruePositiveMin\tTruePositiveMax\tTruePositiveStd\n";
  for (std::vector<OperatingPoint>::const_iterator itr = vCurve.begin(); itr != vCurve.end(); ++itr) {
    os << itr->fThreshold << '\t' << itr->fTruePositiveRate << '\t' << itr->fFalsePositiveCount << '\t' << 
      itr->fThresholdMin << '\t' << itr->fThresholdMax << '\t' << itr->fThresholdStd << '\t' << itr->fTruePositiveRateMin << '\t' << itr->fTruePositiveRateMax << '\t' <<
      itr->fTruePositiveRateStd << '\n';
  }
}

std::vector<OperatingPoint> ResampleCurve(const std::vector<OperatingPoint> &vCurve, const std::set<float> &sFPCounts) {
  typedef PiecewiseLinearInterpolation<float> InterpolatorType;

  InterpolatorType clTPFunc, clTHFunc;

  for (size_t i = 0; i < vCurve.size(); ++i) {
    float x = vCurve[i].fFalsePositiveCount;
    float y = vCurve[i].fTruePositiveRate;

    bool bReplaced = false; // Indicate whether the threshold needs to be replaced

    std::pair<InterpolatorType::Iterator,bool> clResult = clTPFunc.Insert(x, y);
    if (!clResult.second && clResult.first->y < y) {
      // Take the largest true positive rate to be fair (the threshold will be replaced too)
      clTPFunc.Erase(clResult.first);
      clTPFunc.Insert(x, y);
      bReplaced = true;
    }

    y = vCurve[i].fThreshold;
    clResult = clTHFunc.Insert(x, y);

    if (!clResult.second && bReplaced) {
      clTHFunc.Erase(clResult.first);
      clTHFunc.Insert(x, y);
    }
  }

  std::vector<OperatingPoint> vNewCurve;
  vNewCurve.reserve(sFPCounts.size());

  for (std::set<float>::const_iterator itr = sFPCounts.begin(); itr != sFPCounts.end(); ++itr) {
    const float fFalsePositiveRate = *itr;
    const float fThreshold = clTHFunc(fFalsePositiveRate);
    const float fTruePositiveRate = clTPFunc(fFalsePositiveRate);

    vNewCurve.push_back(OperatingPoint(fThreshold, fTruePositiveRate, fFalsePositiveRate));
  }

  return vNewCurve;
}

std::vector<OperatingPoint> AverageFROCCurves(const std::vector<std::vector<OperatingPoint> > &vCurves) {
  std::set<float> sFPCounts;

  for (size_t i = 0; i < vCurves.size(); ++i) {
    for (size_t j = 0; j < vCurves[i].size(); ++j) {
      sFPCounts.insert(vCurves[i][j].fFalsePositiveCount);
    }
  }

  if (sFPCounts.empty())
    return std::vector<OperatingPoint>();

  // Always insert minimum FP rate (threshold >= 1)
  sFPCounts.insert(0.0f);

  std::vector<std::vector<OperatingPoint> > vResampledCurves(vCurves.size());
  for (size_t i = 0; i < vCurves.size(); ++i)
    vResampledCurves[i] = ResampleCurve(vCurves[i], sFPCounts);

  std::vector<OperatingPoint> vMeanCurve(sFPCounts.size());

  // Now compute statistics
  for (size_t i = 0; i < vResampledCurves.size(); ++i) {
    for (size_t j = 0; j < vResampledCurves[i].size(); ++j) {
      const float fThreshold = vResampledCurves[i][j].fThreshold;
      const float fTruePositiveRate = vResampledCurves[i][j].fTruePositiveRate;
      const float fFalsePositiveCount = vResampledCurves[i][j].fFalsePositiveCount;

      OperatingPoint &stMeanPoint = vMeanCurve[j];

      stMeanPoint.fFalsePositiveCount = fFalsePositiveCount;

      const float fThresholdDelta = fThreshold - stMeanPoint.fThreshold;
      const float fTruePositiveRateDelta = fTruePositiveRate - stMeanPoint.fTruePositiveRate;

      stMeanPoint.fThreshold += fThresholdDelta/(i + 1);
      stMeanPoint.fThresholdStd += fThresholdDelta*(fThreshold - stMeanPoint.fThreshold);
  
      stMeanPoint.fTruePositiveRate += fTruePositiveRateDelta/(i + 1);
      stMeanPoint.fTruePositiveRateStd += fTruePositiveRateDelta*(fTruePositiveRate - stMeanPoint.fTruePositiveRate); 

      if (i == 0) {
        stMeanPoint.fThresholdMin = stMeanPoint.fThresholdMax = fThreshold;
        stMeanPoint.fTruePositiveRateMin = stMeanPoint.fTruePositiveRateMax = fTruePositiveRate;
      }
      else {
        stMeanPoint.fThresholdMin = std::min(stMeanPoint.fThresholdMin, fThreshold);
        stMeanPoint.fThresholdMax = std::max(stMeanPoint.fThresholdMax, fThreshold);

        stMeanPoint.fTruePositiveRateMin = std::min(stMeanPoint.fTruePositiveRateMin, fTruePositiveRate);
        stMeanPoint.fTruePositiveRateMax = std::max(stMeanPoint.fTruePositiveRateMax, fTruePositiveRate);
      }
    }
  }

  for (size_t i = 0; i < vMeanCurve.size(); ++i) {
    OperatingPoint &stMeanPoint = vMeanCurve[i];

    if (vResampledCurves.size() < 2) {
      stMeanPoint.fThresholdStd = std::numeric_limits<float>::quiet_NaN();
      stMeanPoint.fTruePositiveRateStd = std::numeric_limits<float>::quiet_NaN();
    }
    else {
      stMeanPoint.fThresholdStd /= vResampledCurves.size()-1;
      stMeanPoint.fTruePositiveRateStd /= vResampledCurves.size()-1;

      stMeanPoint.fThresholdStd = std::sqrt(std::max(0.0f, stMeanPoint.fThresholdStd));
      stMeanPoint.fTruePositiveRateStd = std::sqrt(std::max(0.0f, stMeanPoint.fTruePositiveRateStd));
    }
  }

  return vMeanCurve;
}

OperatingPoint InterpolateOverTPR(const std::vector<OperatingPoint> &vCurve, float fTPR) {
  typedef PiecewiseLinearInterpolation<float> InterpolatorType;

  InterpolatorType clFPFunc;
  InterpolatorType clTHFunc;

  for (size_t i = 0; i < vCurve.size(); ++i) {
    const float x = vCurve[i].fTruePositiveRate;
    float y = vCurve[i].fFalsePositiveCount;

    bool bReplaced = false;

    std::pair<InterpolatorType::Iterator, bool> clResult = clFPFunc.Insert(x, y);
    if (!clResult.second && clResult.first->y > y) {
      clFPFunc.Erase(clResult.first);
      clFPFunc.Insert(x, y);
      bReplaced = true;
    }

    y = vCurve[i].fThreshold;

    clResult = clTHFunc.Insert(x, y);
    if (!clResult.second && bReplaced) {
      clTHFunc.Erase(clResult.first);
      clTHFunc.Insert(x, y);
    }
  }

  const float fThreshold = clTHFunc(fTPR);
  const float fFalsePositiveRate = clFPFunc(fTPR);

  return OperatingPoint(fThreshold, fTPR, fFalsePositiveRate);
}

OperatingPoint InterpolateOverFPC(const std::vector<OperatingPoint> &vCurve, float fFPC) {
  typedef PiecewiseLinearInterpolation<float> InterpolatorType;

  InterpolatorType clTPFunc;
  InterpolatorType clTHFunc;

  for (size_t i = 0; i < vCurve.size(); ++i) {
    const float x = vCurve[i].fFalsePositiveCount;
    float y = vCurve[i].fTruePositiveRate;

    bool bReplaced = false;

    std::pair<InterpolatorType::Iterator, bool> clResult = clTPFunc.Insert(x, y);
    if (!clResult.second && clResult.first->y < y) {
      clTPFunc.Erase(clResult.first);
      clTPFunc.Insert(x, y);
      bReplaced = true;
    }

    y = vCurve[i].fThreshold;

    clResult = clTHFunc.Insert(x, y);
    if (!clResult.second && bReplaced) {
      clTHFunc.Erase(clResult.first);
      clTHFunc.Insert(x, y);
    }
  }

  const float fThreshold = clTHFunc(fFPC);
  const float fTruePositiveRate = clTPFunc(fFPC);

  return OperatingPoint(fThreshold, fTruePositiveRate, fFPC);
}

OperatingPoint InterpolateOverProb(const std::vector<OperatingPoint> &vCurve, float fProb) {
  typedef PiecewiseLinearInterpolation<float> InterpolatorType;

  InterpolatorType clTPFunc;
  InterpolatorType clFPFunc;

  for (size_t i = 0; i < vCurve.size(); ++i) {
    const float x = vCurve[i].fThreshold;
    float y = vCurve[i].fTruePositiveRate;

    bool bReplaced = false;

    std::pair<InterpolatorType::Iterator, bool> clResult = clTPFunc.Insert(x, y);
    if (!clResult.second && clResult.first->y < y) {
      clTPFunc.Erase(clResult.first);
      clTPFunc.Insert(x, y);
      bReplaced = true;
    }

    y = vCurve[i].fFalsePositiveCount;

    clResult = clFPFunc.Insert(x, y);
    if (!clResult.second && bReplaced) {
      clFPFunc.Erase(clResult.first);
      clFPFunc.Insert(x, y);
    }
  }

  const float fFalsePositiveCount = clFPFunc(fProb);
  const float fTruePositiveRate = clTPFunc(fProb);

  return OperatingPoint(fProb, fTruePositiveRate, fFalsePositiveCount);
}