#pragma once

#include <algorithm>
#include <utility>
#include <set>

template<typename RealType>
class PiecewiseLinearInterpolation {
public:
  struct Point {
    RealType x;
    RealType y;

    Point() {
      x = RealType();
      y = RealType();
    }

    Point(RealType x_, RealType y_) {
      x = x_;
      y = y_;
    }

    bool operator<(const Point &stOther) const {
      return x < stOther.x;
    }

    bool operator==(const Point &stOther) const {
      return x == stOther.x;
    }

    bool operator!=(const Point &stOther) const {
      return !(*this == stOther);
    }
  };

  typedef Point PointType;
  typedef std::set<PointType> PointSetType;
  typedef typename PointSetType::iterator Iterator;
  typedef typename PointSetType::const_iterator ConstIterator;
  typedef typename PointSetType::reverse_iterator ReverseIterator;
  typedef typename PointSetType::const_reverse_iterator ConstReverseIterator;

  PiecewiseLinearInterpolation() { }

  template<typename PointIteratorType>
  PiecewiseLinearInterpolation(PointIteratorType begin, PointIteratorType end) {
    Insert(begin, end);
  }

  template<typename XIteratorType, typename YIteratorType>
  PiecewiseLinearInterpolation(XIteratorType xBegin, XIteratorType xEnd, YIteratorType yBegin) {
    Insert(xBegin, xEnd, yBegin);
  }

  size_t Size() const {
    return m_sPoints.size();
  }

  bool Empty() const {
    return m_sPoints.empty();
  }

  Iterator Begin() {
    return m_sPoints.begin();
  }

  Iterator End() {
    return m_sPoints.end();
  }

  ConstIterator Begin() const {
    return m_sPoints.begin();
  }

  ConstIterator End() const {
    return m_sPoints.end();
  }

  ReverseIterator RBegin() {
    return m_sPoints.rbegin();
  }

  ReverseIterator REnd() {
    return m_sPoints.rend();
  }

  ConstReverseIterator RBegin() const {
    return m_sPoints.rbegin();
  }

  ConstReverseIterator REnd() const {
    return m_sPoints.rend();
  }

  std::pair<Iterator,bool> Insert(const PointType &stOther) {
    return m_sPoints.insert(stOther);
  }

  std::pair<Iterator,bool> Insert(const RealType &x, const RealType &y) {
    return Insert(PointType(x,y));
  }

  template<typename PointIteratorType>
  void Insert(PointIteratorType begin, PointIteratorType end) {
    m_sPoints.insert(begin, end);
  }

  template<typename XIteratorType, typename YIteratorType>
  void Insert(XIteratorType xBegin, XIteratorType xEnd, YIteratorType yBegin) {
    for ( ; xBegin != xEnd; ++xBegin, ++yBegin)
      Insert(RealType(*xBegin), RealType(*yBegin));
  }

  void Erase(Iterator itr) {
    m_sPoints.erase(itr);
  }

  void Clear() {
    m_sPoints.clear();
  }

  std::pair<RealType, RealType> Domain() const {
    if (Empty())
      return std::make_pair(RealType(), RealType());

    return std::make_pair(Begin()->x, RBegin()->x);
  }

  std::pair<RealType, RealType> Range() const {
    if (Empty())
      return std::make_pair(RealType(), RealType());

    ConstIterator itr = Begin();

    std::pair<RealType, RealType> stRange(itr->y, itr->y);

    ++itr;

    for ( ; itr != End(); ++itr) {
      stRange.first = std::min(stRange.first, itr->y);
      stRange.second = std::max(stRange.second, itr->y);
    }

    return stRange;
  }

  const PointSetType & Points() const {
    return m_sPoints;
  }

  RealType operator()(const RealType &x) const {
    if (Empty())
      return RealType();

    if (x <= Begin()->x)
      return Begin()->y;

    if (x >= RBegin()->x)
      return RBegin()->y;

    ConstIterator itrPrev = Begin();
    ConstIterator itr = itrPrev;
    ++itr;

    while (itr != End()) {
      if (itrPrev->x <= x && x < itr->x) {
        const RealType t = (itr->x - x)/(itr->x - itrPrev->x);
        return (RealType(1) - t)*(itr->y) + t*(itrPrev->y);
      }

      itrPrev = itr;
      ++itr;
    }

    return RealType(); // Not reached
  }

  RealType Area() const {
    if (Empty())
      return RealType();

    RealType sum = RealType(0);

    ConstIterator itrPrev = Begin();
    ConstIterator itr = itrPrev;
    ++itr;

    while (itr != End()) {
      sum += (itr->x - itrPrev->x)*(itr->y + itrPrev->y);

      itrPrev = itr;
      ++itr;
    }

    return sum/RealType(2);
  }

private:
  PointSetType m_sPoints;
};