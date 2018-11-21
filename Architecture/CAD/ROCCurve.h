#ifndef ROCCURVE_H
#define ROCCURVE_H

#include <iostream>
#include <algorithm>

namespace nih {

class ROCCurve {
public:
  enum { NUM_BINS = 1000 };

  ROCCurve() {
    SetMinMax(0.0f, 1.0f);
    Reset();
  }

  void SetMinMax(float fMinThreshold, float fMaxThreshold) {
    m_fMinThreshold = fMinThreshold;
    m_fMaxThreshold = fMaxThreshold;
  }

  void Reset() {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        std::fill(m_a_counts[i][j], m_a_counts[i][j] + NUM_BINS, (size_t)0);
  }

  float GetMin() const {
    return m_fMinThreshold;
  }

  float GetMax() const {
    return m_fMaxThreshold;
  }

  unsigned int GetNumBins() const {
    return (unsigned int)NUM_BINS;
  }

  float GetDelta() const {
    return (m_fMaxThreshold - m_fMinThreshold)/GetNumBins();
  }

  float GetThresholdByIndex(unsigned int uiIndex) const {
    return GetMin() + GetDelta() * uiIndex;
  }

  float ComputeAUC() const;

  // These give 
  unsigned int GetIndexByThreshold(float fThreshold) const {
    const int iIndex = (int)((fThreshold - m_fMinThreshold)/GetDelta());
    return std::max(0, std::min(iIndex, (int)GetNumBins()-1));
  }

  unsigned int GetIndexByTruePositiveRate(float fTruePositiveRate) const;
  unsigned int GetIndexByFalsePositiveRate(float fFalsePositiveRate) const;
  unsigned int GetIndexByPrecision(float fPrecision) const;
  unsigned int GetRecallByIndex(float fRecall) const {
    return GetIndexByTruePositiveRate(fRecall);
  }

  unsigned int GetIndexBySensitivity(float fSensitivity) const {
    return GetIndexByTruePositiveRate(fSensitivity);
  }

  unsigned int GetIndexBySpecificity(float fSpecificity) const {
    return GetIndexByFalsePositiveRate(fSpecificity);
  }

  size_t GetTruePositiveCountByIndex(unsigned int uiIndex) const {
    return uiIndex < GetNumBins() ? m_a_counts[1][1][uiIndex] : 0;
  }

  size_t GetFalsePositiveCountByIndex(unsigned int uiIndex) const {
    return uiIndex < GetNumBins() ? m_a_counts[0][1][uiIndex] : 0;
  }

  size_t GetTrueNegativeCountByIndex(unsigned int uiIndex) const {
    return uiIndex < GetNumBins() ? m_a_counts[0][0][uiIndex] : GetNumNegatives();
  }

  size_t GetFalseNegativeCountByIndex(unsigned int uiIndex) const {
    return uiIndex < GetNumBins() ? m_a_counts[1][0][uiIndex] : GetNumPositives();
  }

  size_t GetNumPositives() const {
    return GetTruePositiveCountByIndex(0) + GetFalseNegativeCountByIndex(0);
  }

  size_t GetNumNegatives() const {
    return GetTrueNegativeCountByIndex(0) + GetFalsePositiveCountByIndex(0);
  }

  float GetMaxPrecision() const {
    const size_t posPlusNeg = GetNumPositives() + GetNumNegatives();
    return posPlusNeg > 0 ? GetNumPositives()/(float)posPlusNeg : 0.0f;
  }

  float GetSensitivityByIndex(unsigned int uiIndex) const {
    return GetTruePositiveRateByIndex(uiIndex);
  }

  float GetSpecificityByIndex(unsigned int uiIndex) const {
    return GetFalsePositiveRateByIndex(uiIndex);
  }

  float GetPrecisionByIndex(unsigned int uiIndex) const {
    const size_t tpPlusFp = GetTruePositiveCountByIndex(uiIndex) + GetFalsePositiveCountByIndex(uiIndex);
    return tpPlusFp > 0 ? GetTruePositiveCountByIndex(uiIndex)/(float)tpPlusFp : 0.0f;
  }

  float GetRecallByIndex(unsigned int uiIndex) const {
    return GetTruePositiveRateByIndex(uiIndex);
  }

  float GetTruePositiveRateByIndex(unsigned int iIndex) const {
    return GetNumPositives() > 0 ? GetTruePositiveCountByIndex(iIndex)/(float)GetNumPositives() : 0.0f;
  }

  float GetFalsePositiveRateByIndex(unsigned int uiIndex) const {
    return GetNumNegatives() > 0 ? GetFalsePositiveCountByIndex(uiIndex)/(float)GetNumNegatives() : 0.0f;
  }
  
  // These do linear interpolation
  float GetTruePositiveRateByThreshold(float fThreshold) const;
  float GetFalsePositiveRateByThreshold(float fThreshold) const;

  float GetSensitivityByThreshold(float fThreshold) const {
    return GetTruePositiveRateByThreshold(fThreshold);
  }

  float GetSpecificityByThreshold(float fThreshold) const {
    return GetFalsePositiveRateByThreshold(fThreshold);
  }

  float GetPrecisionByThreshold(float fThreshold) const;

  float GetRecallByThreshold(float fThreshold) const {
    return GetTruePositiveRateByThreshold(fThreshold);
  }

  float GetThresholdByTruePositiveRate(float fTruePositiveRate) const;
  float GetThresholdByFalsePositiveRate(float fFalsePositiveRate) const;

  float GetThresholdBySensitivity(float fSensitivity) const {
    return GetThresholdByTruePositiveRate(fSensitivity);
  }

  float GetThresholdBySpecificity(float fSpecificity) const {
    return GetThresholdByFalsePositiveRate(fSpecificity);
  }

  float GetThresholdByPrecision(float fPrecision) const;

  float GetThresholdByRecall(float fRecall) const {
    return GetThresholdByFalsePositiveRate(fRecall);
  }

  void Update(float fProbability, int iGTLabel) {
    iGTLabel = (iGTLabel > 0);

    for (unsigned int i = 0; i < GetNumBins(); ++i) {
      const float fThreshold = GetThresholdByIndex(i);
      //const int iPredictedLabel = (fProbability > fThreshold);
      const int iPredictedLabel = (fProbability > fThreshold);

      ++m_a_counts[iGTLabel][iPredictedLabel][i];
    }
  }

  void MergeCounts(const ROCCurve &clOther) {
    for (unsigned int i = 0; i < GetNumBins(); ++i) {
      m_a_counts[0][0][i] += clOther.m_a_counts[0][0][i];
      m_a_counts[0][1][i] += clOther.m_a_counts[0][1][i];
      m_a_counts[1][0][i] += clOther.m_a_counts[1][0][i];
      m_a_counts[1][1][i] += clOther.m_a_counts[1][1][i];
    }
  }

private:
  float m_fMinThreshold, m_fMaxThreshold;
  size_t m_a_counts[2][2][NUM_BINS]; // Actual x Predicted x Bins
};

void PrintROC(std::ostream &os, const ROCCurve &clROC);
void PrintPrecisionRecall(std::ostream &os, const ROCCurve &clROC);
void PrintConfusionMatrix(std::ostream &os, const ROCCurve &clROC);

class ROCCurveFormat {
public:
  const ROCCurve &clROC;

  explicit ROCCurveFormat(const ROCCurve &clROC_)
  : clROC(clROC_) { }

  virtual ~ROCCurveFormat() { }

  virtual void Print(std::ostream &os) const = 0;
};

template<void (*PrintFunction)(std::ostream &, const ROCCurve &)>
class ROCCurveFormatTemplate : public ROCCurveFormat {
public:
  explicit ROCCurveFormatTemplate(const ROCCurve &clROC_)
  : ROCCurveFormat(clROC_) { }

  void Print(std::ostream &os) const {
    (*PrintFunction)(os, clROC);
  }
};

class FROCCurveFormat : public ROCCurveFormat {
public:
  const unsigned int uiNumCases;

  explicit FROCCurveFormat(const ROCCurve &clROC_, unsigned int uiNumCases_ = 1)
  : ROCCurveFormat(clROC_), uiNumCases(uiNumCases_) { }

  void Print(std::ostream &os) const;
};

typedef ROCCurveFormatTemplate<&PrintROC> ROCFormat;
typedef ROCCurveFormatTemplate<&PrintPrecisionRecall> PrecisionRecallFormat;
typedef ROCCurveFormatTemplate<&PrintConfusionMatrix> ConfusionMatrixFormat;

inline std::ostream & operator<<(std::ostream &os, const ROCCurveFormat &clROCFormat) {
  clROCFormat.Print(os);
  return os;
}

inline std::ostream & operator<<(std::ostream &os, const ROCCurve &clROC) {
  return os << nih::ROCFormat(clROC);
}

} // end namespace nih



#endif // ROCCURVE_H
