#include "ROCCurve.h"

namespace nih {

float ROCCurve::ComputeAUC() const {
  float fSum = 0.0f;

  // Threshold 0: TPR = 1, FPR = 1
  fSum += (1.0f - GetFalsePositiveRateByIndex(0)) * (1.0f + GetTruePositiveRateByIndex(0));

  for (unsigned int i = 0; i < GetNumBins()-1; ++i) {
    const float fDeltaX = GetFalsePositiveRateByIndex(i) - GetFalsePositiveRateByIndex(i+1);
    const float fDeltaY = GetTruePositiveRateByIndex(i) + GetTruePositiveRateByIndex(i+1);
    fSum += fDeltaX*fDeltaY;
  }

  // Threshold 1: TPR = 0, FPR = 0
  fSum += (GetFalsePositiveRateByIndex(GetNumBins()-1) - 0.0f) * (0.0f + GetTruePositiveRateByIndex(GetNumBins()-1));

  return 0.5f*fSum;
}

unsigned int ROCCurve::GetIndexByTruePositiveRate(float fTruePositiveRate) const {
  if (fTruePositiveRate >= 1.0f)
    return 0;

  if (fTruePositiveRate <= 0.0f)
    return GetNumBins();

  for (unsigned int i = 1; i < GetNumBins(); ++i) {
    if (GetTruePositiveRateByIndex(i) <= fTruePositiveRate)
      return i-1;
  }

  return 0; // Uh?
}

unsigned int ROCCurve::GetIndexByFalsePositiveRate(float fFalsePositiveRate) const {
  if (fFalsePositiveRate >= 1.0f)
    return 0;

  if (fFalsePositiveRate <= 0.0f)
    return GetNumBins();

  for (unsigned int i = 1; i < GetNumBins(); ++i) {
    if (GetFalsePositiveRateByIndex(i) <= fFalsePositiveRate)
      return i-1;
  }

  return 0; // Uh? 
}

unsigned int ROCCurve::GetIndexByPrecision(float fPrecision) const {
  if (fPrecision >= GetMaxPrecision())
    return 0;

  if (fPrecision <= 0.0f)
    return GetNumBins();

  for (unsigned int i = 1; i < GetNumBins(); ++i) {
    if (GetPrecisionByIndex(i) <= fPrecision)
      return i-1;
  }

  return 0; // Uh?
}

float ROCCurve::GetTruePositiveRateByThreshold(float fThreshold) const {
  if (fThreshold < GetMin())
    return 1.0f;

  if (fThreshold >= GetMax())
    return 0.0f;

  const unsigned int uiIndex = GetIndexByThreshold(fThreshold);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fTruePositiveRate1 = GetTruePositiveRateByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fTruePositiveRate2 = GetTruePositiveRateByIndex(uiIndex+1);

  const float t = (fThreshold - fThreshold1)/(fThreshold2 - fThreshold1);

  return fTruePositiveRate1 + (fTruePositiveRate2 - fTruePositiveRate1) * t;
}

float ROCCurve::GetFalsePositiveRateByThreshold(float fThreshold) const {
  if (fThreshold < GetMin())
    return 1.0f;

  if (fThreshold >= GetMax())
    return 0.0f;

  const unsigned int uiIndex = GetIndexByThreshold(fThreshold);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fFalsePositiveRate1 = GetFalsePositiveRateByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fFalsePositiveRate2 = GetFalsePositiveRateByIndex(uiIndex+1);

  const float t = (fThreshold - fThreshold1)/(fThreshold2 - fThreshold1);

  return fFalsePositiveRate1 + (fFalsePositiveRate2 - fFalsePositiveRate1) * t;
}

float ROCCurve::GetPrecisionByThreshold(float fThreshold) const {
  if (fThreshold < GetMin())
    return GetMaxPrecision();

  if (fThreshold >= GetMax())
    return 0.0f;

  const unsigned int uiIndex = GetIndexByThreshold(fThreshold);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fPrecision1 = GetPrecisionByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fPrecision2 = GetPrecisionByIndex(uiIndex+1);

  const float t = (fThreshold - fThreshold1)/(fThreshold2 - fThreshold1);

  return fPrecision1 + (fPrecision2 - fPrecision1) * t;
}

float ROCCurve::GetThresholdByTruePositiveRate(float fTruePositiveRate) const {
  if (fTruePositiveRate >= 1.0f)
    return GetMin();

  if (fTruePositiveRate <= 0.0f)
    return GetMax();

  const unsigned int uiIndex = GetIndexByTruePositiveRate(fTruePositiveRate);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fTruePositiveRate1 = GetTruePositiveRateByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fTruePositiveRate2 = GetTruePositiveRateByIndex(uiIndex+1);

  const float t = (fTruePositiveRate - fTruePositiveRate1)/(fTruePositiveRate2 - fTruePositiveRate1);

  return fThreshold1 + (fThreshold2 - fThreshold1) * t;
}

float ROCCurve::GetThresholdByFalsePositiveRate(float fFalsePositiveRate) const {
  if (fFalsePositiveRate >= 1.0f)
    return GetMin();

  if (fFalsePositiveRate <= 0.0f)
    return GetMax();

  const unsigned int uiIndex = GetIndexByFalsePositiveRate(fFalsePositiveRate);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fFalsePositiveRate1 = GetFalsePositiveRateByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fFalsePositiveRate2 = GetFalsePositiveRateByIndex(uiIndex+1);

  const float t = (fFalsePositiveRate - fFalsePositiveRate1)/(fFalsePositiveRate2 - fFalsePositiveRate1);

  return fThreshold1 + (fThreshold2 - fThreshold1) * t;
}

float ROCCurve::GetThresholdByPrecision(float fPrecision) const {
  if (fPrecision >= GetMaxPrecision())
    return GetMin();

  if (fPrecision <= 0.0f)
    return GetMax();

  const unsigned int uiIndex = GetIndexByPrecision(fPrecision);

  const float fThreshold1 = GetThresholdByIndex(uiIndex);
  const float fPrecision1 = GetPrecisionByIndex(uiIndex);
  const float fThreshold2 = GetThresholdByIndex(uiIndex+1);
  const float fPrecision2 = GetPrecisionByIndex(uiIndex+1);

  const float t = (fPrecision - fPrecision1)/(fPrecision2 - fPrecision1);

  return fThreshold1 + (fThreshold2 - fThreshold1) * t;
}

void PrintROC(std::ostream &os, const nih::ROCCurve &clROC) {
  os << "# Threshold\tTrue Positive Rate\tFalse Positive Rate\n";

  for (unsigned int i = 0; i < clROC.GetNumBins(); ++i) {
    const float fThreshold = clROC.GetThresholdByIndex(i);
    const float fTruePositiveRate = clROC.GetTruePositiveRateByIndex(i);
    const float fFalsePositiveRate = clROC.GetFalsePositiveRateByIndex(i);

    os << fThreshold << '\t' << fTruePositiveRate << '\t' << fFalsePositiveRate << '\n';
  }

}

void PrintPrecisionRecall(std::ostream &os, const nih::ROCCurve &clROC) {

  os << "# Threshold\tPrecision\tRecall\n";

  for (unsigned int i = 0; i < clROC.GetNumBins(); ++i) {
    const float fThreshold = clROC.GetThresholdByIndex(i);
    const float fPrecision = clROC.GetPrecisionByIndex(i);
    const float fRecall = clROC.GetRecallByIndex(i);

    os << fThreshold << '\t' << fPrecision << '\t' << fRecall << '\n';
  }
}

void PrintConfusionMatrix(std::ostream &os, const nih::ROCCurve &clROC) {

  os << "# Threshold\tTrue Positives\tFalse Positives\tTrue Negatives\tFalse Negatives\n";

  for (unsigned int i = 0; i < clROC.GetNumBins(); ++i) {
    const float fThreshold = clROC.GetThresholdByIndex(i);
    const size_t tpCount = clROC.GetTruePositiveCountByIndex(i);
    const size_t fpCount = clROC.GetFalsePositiveCountByIndex(i);
    const size_t tnCount = clROC.GetTrueNegativeCountByIndex(i);
    const size_t fnCount = clROC.GetFalseNegativeCountByIndex(i);

    os << fThreshold << '\t' << tpCount << '\t' << fpCount << '\t' << tnCount << '\t' << fnCount << '\n';
  }
}

void FROCCurveFormat::Print(std::ostream &os) const {
  os << "# Threshold\tTrue Positive Rate\tAverage False Positive Count\n";

  for (unsigned int i = 0; i < clROC.GetNumBins(); ++i) {
    const float fThreshold = clROC.GetThresholdByIndex(i);
    const float fTruePositiveRate = clROC.GetTruePositiveRateByIndex(i);
    const size_t falsePositiveCount = clROC.GetFalsePositiveCountByIndex(i);

    os << fThreshold << '\t' << fTruePositiveRate << '\t' << (falsePositiveCount / (float)uiNumCases) << '\n';
  }
}

} // end namespace nih
