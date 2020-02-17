/*
 * Example Code for the Paper
 *
 * "Mesh Saliency via Local Curvature Entropy", Eurographics 2016
 *
 * by M. Limper, A. Kuijper and D. Fellner
 *
 */

#ifndef LCE_FLOAT_HISTOGRAM_H
#define LCE_FLOAT_HISTOGRAM_H


#include <iostream>
#include <map>


namespace lce
{

//------------------------------------------------------------------------------------------------------------------

class FloatHistogram
{
public:

    FloatHistogram();


    FloatHistogram(const FloatHistogram & other);


    void add(float value, float weight);


    float getPercentile(float p) const;


    float getMin() const
    {
        return _min;
    }


    float getMax() const
    {
        return _max;
    }


    const std::map<float, float> & getValues() const
    {
        return _count;
    }


private:

    std::map<float, float> _count;

    float _totalCount;

    float _min;

    float _max;


}; //class FloatHistogram

//------------------------------------------------------------------------------------------------------------------

} //namespace lce


#endif

