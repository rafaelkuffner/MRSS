/*
 * Example Code for the Paper
 *
 * "Mesh Saliency via Local Curvature Entropy", Eurographics 2016
 *
 * by M. Limper, A. Kuijper and D. Fellner
 *
 */

#include "Histogram.h"

#include <float.h>
#include <algorithm>


namespace lce
{

//------------------------------------------------------------------------------------------------------------------

Histogram::Histogram()
    : _totalCount(0), _min(FLT_MAX), _max(-FLT_MAX)
{
}

//------------------------------------------------------------------------------------------------------------------

void Histogram::add(float value)
{
    _min = std::min(value, _min);
    _max = std::max(value, _max);

    ++_count[value];
    ++_totalCount;
}

//------------------------------------------------------------------------------------------------------------------

float Histogram::getPercentile(float p) const
{
    if (p == 0.0f)
        return _min;

    if (p == 1.0f)
        return _max;

	unsigned n = 0;
	for (auto mIt = _count.begin(); mIt != _count.end(); ++mIt)
	{
		if (static_cast<float>(n) / _totalCount >= p)
		{
			return mIt->first;
		}
		n += mIt->second;
	}

    return _max;
}

float Histogram::getMinNoOutliers() const {
	float median = getPercentile(0.5);
	float iqr = getPercentile(0.75) - getPercentile(0.25);
	float minCapped = median - (iqr * 15);
	minCapped = getMin() < minCapped ? minCapped : getMin();

	return minCapped;

}

float Histogram::getMaxNoOutliers() const {
	float median = getPercentile(0.5);
	float iqr = getPercentile(0.75) - getPercentile(0.25);
	float maxCapped = median + (iqr * 15);
	maxCapped = getMax() > maxCapped ? maxCapped : getMax();

	return maxCapped;
}
//------------------------------------------------------------------------------------------------------------------

}
