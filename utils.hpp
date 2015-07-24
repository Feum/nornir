/*
 * utils.hpp
 *
 * Created on: 09/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of AdaptiveFastFlow.
 *
 *  AdaptiveFastFlow is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  AdaptiveFastFlow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with AdaptiveFastFlow.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

/*!
 * \file utils.hpp
 * \brief Implementation of various utilities.
 **/

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>

namespace adpff{

/**
 * Represents a moving average technique.
 * Requirement: There must exists a method 'T squareRoot(const T&)'.
 */
template <typename T> class MovingAverage{
public:
    virtual ~MovingAverage(){;}

    virtual void add(const T& value) = 0;

    virtual T getLastSample() const = 0;

    virtual void reset() = 0;

    virtual size_t size() const = 0;

    virtual T average() const = 0;

    virtual T variance() const = 0;

    virtual T standardDeviation() const = 0;

    virtual T coefficientVariation() const = 0;
};

template<typename T> class MovingAverageSimple: public MovingAverage<T>{
    template<typename V>
    friend std::ostream& operator<<(std::ostream& os,
                                    const MovingAverageSimple<V>& obj);
private:
    std::vector<T> _windowImpl;
    size_t _span;
    size_t _nextIndex;
    size_t _storedValues;
    T _lastSample;
    T _oldAverage, _average;
    T _tmpVariance;
    T _oldVariance, _variance;
    T _standardDeviation;
    T _coefficientVariation;
public:
    MovingAverageSimple(size_t span):_span(span), _nextIndex(0),
                                     _storedValues(0){
        _windowImpl.resize(_span);
    }

    void add(const T& value){
        _lastSample = value;
        if(_storedValues == 0){
            ++_storedValues;
            _oldAverage = _average = value;
        }else if(_storedValues < _span){
            ++_storedValues;
            _average = _oldAverage + (value - _oldAverage)/_storedValues;
            _tmpVariance = _tmpVariance + (value - _oldAverage)*
                                          (value - _average);
            _variance = _tmpVariance / _storedValues;

            _oldAverage = _average;
            _oldVariance = _variance;
        }else{
            T removedValue = _windowImpl.at(_nextIndex);
            _average = _oldAverage + (value - _windowImpl.at(_nextIndex))/
                                         _span;
            _variance = _oldVariance + (value - _average +
                                           removedValue - _oldAverage)*
                                          (value - removedValue)/(_span-1);

            _windowImpl.at(_nextIndex) = value;
            _nextIndex = (_nextIndex + 1) % _span;
        }
        _standardDeviation = squareRoot(_variance);
        _coefficientVariation = (_standardDeviation / _average) * 100.0;
    }

    T getLastSample() const{
        return _lastSample;
    }

    void reset(){
        _windowImpl.clear();
        _windowImpl.resize(_span);
        _nextIndex = 0;
        _storedValues = 0;
    }

    size_t size() const{
        return _storedValues;
    }

    T average() const{
        return _average;
    }

    T variance() const{
        return _variance;
    }

    T standardDeviation() const{
        return _standardDeviation;
    }

    T coefficientVariation() const{
        return _coefficientVariation;
    }
};

template<class T>
std::ostream& operator<<(std::ostream& os, const MovingAverage<T>& obj){
    os << "[";
    os << "Last sample: " << obj.getLastSample();
    os << "Average: " << obj.average();
    os << "Variance: " << obj.variance();
    os << "StdDev: " << obj.standardDeviation();
    os << "CoeffVar: " << obj.coefficientVariation();
    os << "]";
    return os;
}

}

#endif /* UTILS_HPP_ */
