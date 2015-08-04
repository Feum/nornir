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

#define NSECS_IN_SECS 1000000000.0

namespace adpff{

/**
 * Represents a moving average technique.
 * Requirement: There must exists the following functions:
 *   - 'T squareRoot(const T&)' to compute the square root.
 *   - 'void regularize(T&)' to set the values < 0 to zero.
 */
template <typename T> class Smoother{
    template<typename V>
    friend std::ostream& operator<<(std::ostream& os,
                                    const Smoother<V>& obj);

    template<typename V>
    friend std::ofstream& operator<<(std::ofstream& os,
                                     const Smoother<V>& obj);
public:
    virtual ~Smoother(){;}

    virtual double getSmoothingFactor() const = 0;

    virtual void setSmoothingFactor(double s) = 0;

    virtual void add(const T& value) = 0;

    virtual T getLastSample() const = 0;

    virtual void reset() = 0;

    virtual size_t size() const = 0;

    virtual T average() const = 0;

    virtual T variance() const = 0;

    virtual T standardDeviation() const = 0;

    virtual T coefficientVariation() const = 0;
};

template<typename T> class MovingAverageSimple: public Smoother<T>{
private:
    std::vector<T> _windowImpl;
    size_t _span;
    size_t _nextIndex;
    size_t _storedValues;
    T _lastSample;
    T _oldAverage, _average;
    T _oldTmpVariance, _tmpVariance;
    T _oldVariance, _variance;
    T _standardDeviation;
    T _coefficientVariation;
public:
    MovingAverageSimple(size_t span):_span(span), _nextIndex(0),
                                     _storedValues(0){
        _windowImpl.resize(_span);
    }

    double getSmoothingFactor() const{
        return _span;
    }

    void setSmoothingFactor(double s){
        reset();
        _span = s;
        _windowImpl.resize(_span);
    }

    void add(const T& value){
        _lastSample = value;
        if(_storedValues == 0){
            ++_storedValues;
            _oldAverage = _average = value;
        }else if(_storedValues < _span){
            ++_storedValues;
            _average = _oldAverage + (value - _oldAverage) / _storedValues;
            _tmpVariance = _oldTmpVariance + (value - _oldAverage)*
                                             (value - _average);
            _variance = _tmpVariance / (_storedValues - 1);
        }else{
            T removedValue = _windowImpl.at(_nextIndex);
            _average = _oldAverage + (value - removedValue)/
                                      _span;
            _variance = _oldVariance + (value - _average +
                                        removedValue - _oldAverage)*
                                        (value - removedValue)/(_span - 1);
        }

        _oldAverage = _average;
        _oldTmpVariance = _tmpVariance;
        _oldVariance = _variance;

        regularize(_variance);

        _standardDeviation = squareRoot(_variance);
        _coefficientVariation = (_standardDeviation / _average) * 100.0;

        _windowImpl.at(_nextIndex) = value;
        _nextIndex = (_nextIndex + 1) % _span;
    }

    T getLastSample() const{
        return _lastSample;
    }

    void reset(){
        _windowImpl.clear();
        _windowImpl.resize(_span);
        _nextIndex = 0;
        _storedValues = 0;
        _lastSample = T();
        _oldAverage = T();
        _average = T();
        _oldTmpVariance = T();
        _tmpVariance = T();
        _oldVariance = T();
        _variance = T();
        _standardDeviation = T();
        _coefficientVariation = T();
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

template<typename T> class MovingAverageExponential: public Smoother<T>{
private:
    double _alpha;
    size_t _storedValues;
    T _lastSample;
    T _average;
    T _variance;
    T _standardDeviation;
    T _coefficientVariation;
public:
    MovingAverageExponential(double alpha):_alpha(alpha),_storedValues(0){
        if(_alpha < 0 || _alpha > 1.0){
            throw std::runtime_error("Alpha must be between 0 and 1 "
                                     "(included)");
        }
    }

    double getSmoothingFactor() const{
        return _alpha;
    }

    void setSmoothingFactor(double s){
        _alpha = s;
    }

    void add(const T& value){
        ++_storedValues;
        _lastSample = value;
        if(_storedValues == 1){
            _average = value;
        }else{
            T diff = value - _average;
            T incr = diff * _alpha;
            _average += incr;
            _variance = (_variance + diff * incr) * (1 - _alpha);
        }
        regularize(_variance);
        _standardDeviation = squareRoot(_variance);
        _coefficientVariation = (_standardDeviation / _average) * 100.0;
    }

    T getLastSample() const{
        return _lastSample;
    }

    void reset(){
        _storedValues = 0;
        _lastSample = T();
        _average = T();
        _variance = T();
        _standardDeviation = T();
        _coefficientVariation = T();
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
std::ostream& operator<<(std::ostream& os, const Smoother<T>& obj){
    os << "==============================" << std::endl;
    os << "Last sample: " << obj.getLastSample() << std::endl;
    os << "Average: " << obj.average() << std::endl;
    os << "CoeffVar: " << obj.coefficientVariation() << std::endl;
    os << "Variance: " << obj.variance() << std::endl;
    os << "StdDev: " << obj.standardDeviation() << std::endl;
    os << "==============================" << std::endl;
    return os;
}

template<class T>
std::ofstream& operator<<(std::ofstream& os, const Smoother<T>& obj){
    os << obj.getLastSample() << "\t";
    os << obj.average() << "\t";
    os << obj.coefficientVariation() << "\t";
    os << obj.variance() << "\t";
    os << obj.standardDeviation() << "\t";
    return os;
}

}

#endif /* UTILS_HPP_ */
