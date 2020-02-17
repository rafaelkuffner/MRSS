/*
 * Example Code for the Paper
 *
 * "Mesh Saliency via Local Curvature Entropy", Eurographics 2016
 *
 * by M. Limper, A. Kuijper and D. Fellner
 *
 */

#include "FloatHistogram.h"

#include <float.h>
#include <algorithm>


namespace lce
{

//------------------------------------------------------------------------------------------------------------------

FloatHistogram::FloatHistogram()
    : _totalCount(0.0f), _min(FLT_MAX), _max(-FLT_MAX)
{
}

//------------------------------------------------------------------------------------------------------------------

FloatHistogram::FloatHistogram(const FloatHistogram & other)
{
    this->_min  = other._min;
    this->_max  = other._max;
    this->_count      = other._count;
    this->_totalCount = other._totalCount;
}

//------------------------------------------------------------------------------------------------------------------

void FloatHistogram::add(float value, float weight)
{
    _min = std::min(value, _min);
    _max = std::max(value, _max);

    _count[value] += weight;
    _totalCount   += weight;
}

//------------------------------------------------------------------------------------------------------------------

float FloatHistogram::getPercentile(float p) const
{
    if (p == 0.0f)
        return _min;

    if (p == 1.0f)
        return _max;

    std::map<float, float>::const_iterator mIt, mEnd;

    float n = 0;

    for (mIt = _count.begin(), mEnd = _count.end(); mIt != mEnd; ++mIt)
    {
        if (n / _totalCount >= p)
        {
            return mIt->first;
        }

        n += mIt->second;
    }

    return _max;
}

//------------------------------------------------------------------------------------------------------------------

} //namespace lce
