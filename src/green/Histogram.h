/*
 * Example Code for the Paper
 *
 * "Mesh Saliency via Local Curvature Entropy", Eurographics 2016
 *
 * by M. Limper, A. Kuijper and D. Fellner
 *
 */

#ifndef LCE_HISTOGRAM_H
#define LCE_HISTOGRAM_H


#include <iostream>
#include <map>


namespace lce
{

//------------------------------------------------------------------------------------------------------------------

class Histogram
{
public:

    Histogram();


    void add(float value);


    float getPercentile(float p) const;


    float getMin() const
    {
        return _min;
    }


    float getMax() const
    {
        return _max;
    }

	float getMinNoOutliers() const;

	float getMaxNoOutliers() const;

    const std::map<float, unsigned> & getValues() const
    {
        return _count;
    }



private:

    std::map<float, unsigned> _count;

    unsigned _totalCount;

    float _min;

    float _max;


}; //class Histogram

//------------------------------------------------------------------------------------------------------------------

} //namespace lce


#endif

