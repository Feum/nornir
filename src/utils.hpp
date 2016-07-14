/*
 * utils.hpp
 *
 * Created on: 09/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

/*!
 * \file utils.hpp
 * \brief Implementation of various utilities.
 **/

#ifndef NORNIR_UTILS_HPP_
#define NORNIR_UTILS_HPP_

#include "external/Mammut/mammut/mammut.hpp"

#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <functional>
#include <numeric>

#define MSECS_IN_SECS 1000.0 // Milliseconds in 1 second
#define NSECS_IN_SECS 1000000000.0 // Nanoseconds in 1 second

namespace nornir{

inline double ticksToSeconds(double ticks, double ticksPerNs){
    return (ticks/ticksPerNs)/NSECS_IN_SECS;
}

inline double ticksToMilliseconds(double ticks, double ticksPerNs){
    return ticksToSeconds(ticks, ticksPerNs)*1000;
}

inline double average(const std::vector<double>& v){
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

inline double stddev(const std::vector<double>& v, double average){
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   bind2nd(std::minus<double>(), average));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return sqrt(sq_sum / v.size());
}

inline double stddev(const std::vector<double>& v){
    return stddev(v, average(v));
}

typedef struct MonitoredSample{
    mammut::energy::Joules watts; ///< Consumed watts.
    double utilisation; ///< Utilisation [0, 100].
    double bandwidth; ///< Bandwidth of the entire farm (real).
    double latency; ///< Average latency of a worker (nanoseconds).

    MonitoredSample():watts(0), utilisation(0), bandwidth(0), latency(0){;}

    double getMaximumBandwidth(){
        return bandwidth / (utilisation / 100.0);
    }

    void swap(MonitoredSample& x){
        using std::swap;

        swap(watts, x.watts);
        swap(utilisation, x.utilisation);
        swap(bandwidth, x.bandwidth);
        swap(latency, x.latency);
    }

    MonitoredSample& operator=(MonitoredSample rhs){
        swap(rhs);
        return *this;
    }

    MonitoredSample& operator+=(const MonitoredSample& rhs){
        watts += rhs.watts;
        utilisation += rhs.utilisation;
        bandwidth += rhs.bandwidth;
        latency += rhs.latency;
        return *this;
    }

    MonitoredSample& operator-=(const MonitoredSample& rhs){
        watts -= rhs.watts;
        utilisation -= rhs.utilisation;
        bandwidth -= rhs.bandwidth;
        latency -= rhs.latency;
        return *this;
    }

    MonitoredSample& operator*=(const MonitoredSample& rhs){
        watts *= rhs.watts;
        utilisation *= rhs.utilisation;
        bandwidth *= rhs.bandwidth;
        latency *= rhs.latency;
        return *this;
    }

    MonitoredSample& operator/=(const MonitoredSample& rhs){
        watts /= rhs.watts;
        utilisation /= rhs.utilisation;
        bandwidth /= rhs.bandwidth;
        latency /= rhs.latency;
        return *this;
    }

    MonitoredSample operator/=(double x){
        watts /= x;
        utilisation /= x;
        bandwidth /= x;
        latency /= x;
        return *this;
    }

    MonitoredSample operator*=(double x){
        watts *= x;
        utilisation *= x;
        bandwidth *= x;
        latency *= x;
        return *this;
    }
}MonitoredSample;

inline MonitoredSample operator+(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r += rhs;
    return r;
}

inline MonitoredSample operator-(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r -= rhs;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r *= rhs;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r /= rhs;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r /= x;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r *= x;
    return r;
}

inline std::ostream& operator<<(std::ostream& os, const MonitoredSample& obj){
    os << "[";
    os << "Watts: " << obj.watts << " ";
    os << "BandwidthMax: " << obj.utilisation << " ";
    os << "Bandwidth: " << obj.bandwidth << " ";
    os << "Latency: " << obj.latency << " ";
    os << "]";
    return os;
}

inline std::ofstream& operator<<(std::ofstream& os, const MonitoredSample& obj){
    os << obj.watts << "\t";
    os << obj.utilisation << "\t";
    os << obj.bandwidth << "\t";
    os << obj.latency << "\t";
    return os;
}

inline MonitoredSample squareRoot(const MonitoredSample& x){
    MonitoredSample r;
    r.watts = x.watts?sqrt(x.watts):0;
    r.utilisation = x.utilisation?sqrt(x.utilisation):0;
    r.bandwidth = x.bandwidth?sqrt(x.bandwidth):0;
    r.latency = x.latency?sqrt(x.latency):0;
    return r;
}

inline void regularize(MonitoredSample& x){
    if(x.watts < 0){x.watts = 0;}
    if(x.utilisation < 0){x.utilisation = 0;}
    if(x.bandwidth < 0){x.bandwidth = 0;}
    if(x.latency < 0){x.latency = 0;}
}

inline double squareRoot(const double& x){
    return x?sqrt(x):0;
}

inline void regularize(double& x){
    if(x < 0){x = 0;}
}

/**
 * Represents a smoothing technique.
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

    /**
     * Returns the smoothing factor.
     * @return The smoothing factor.
     */
    virtual double getSmoothingFactor() const = 0;

    /**
     * Sets the smoothing factor.
     * @param s The smoothing factor.
     */
    virtual void setSmoothingFactor(double s) = 0;

    /**
     * Adds a sample to the smoother.
     * @param value The sample to be added.
     */
    virtual void add(const T& value) = 0;

    /**
     * Gets the last stored sample.
     * @return The last stored sample.
     */
    virtual T getLastSample() const = 0;

    /**
     * Resets the smoother.
     */
    virtual void reset() = 0;

    /**
     * Returns the number of samples stored by the smoother.
     * @return The number of samples stored by the smoother.
     */
    virtual size_t size() const = 0;

    /**
     * Returns the average of the stored samples.
     * @return The average of the stored samples.
     */
    virtual T average() const = 0;

    /**
     * Returns the variance of the stored samples.
     * @return The variance of the stored samples.
     */
    virtual T variance() const = 0;

    /**
     * Returns the standard deviation of the stored samples.
     * @return The standard deviation of the stored samples.
     */
    virtual T standardDeviation() const = 0;

    /**
     * Returns the coefficient of variation of the stored samples [0, 100].
     * @return The coefficient of variation of the stored samples [0, 100].
     */
    virtual T coefficientVariation() const = 0;
};

/**
 * Smoothing technique: Moving average
 */
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
            _variance = T();
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

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v){
    if(v.empty()){
        out << '[';
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}

/*****************************************************************************************/
/* To disable compiler warnings.                                                         */
/* Code obtained from: https://github.com/facebook/folly/blob/master/folly/Portability.h */
/*****************************************************************************************/

// Generalize warning push/pop.
#if defined(_MSC_VER)
# define PUSH_WARNING __pragma(warning(push))
# define POP_WARNING __pragma(warning(pop))
// Disable the GCC warnings.
# define GCC_DISABLE_WARNING(warningName)
# define MSVC_DISABLE_WARNING(warningNumber) __pragma(warning(disable: warningNumber))
#elif defined(__clang__) || __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
# define PUSH_WARNING _Pragma("GCC diagnostic push")
# define POP_WARNING _Pragma("GCC diagnostic pop")
#define GCC_DISABLE_WARNING_INTERNAL3(warningName) #warningName
#define GCC_DISABLE_WARNING_INTERNAL2(warningName) \
  GCC_DISABLE_WARNING_INTERNAL3(warningName)
#define GCC_DISABLE_WARNING(warningName)                       \
  _Pragma(GCC_DISABLE_WARNING_INTERNAL2(GCC diagnostic ignored \
          GCC_DISABLE_WARNING_INTERNAL3(-W##warningName)))
// Disable the MSVC warnings.
# define MSVC_DISABLE_WARNING(warningNumber)
#else
# define PUSH_WARNING
# define POP_WARNING
# define GCC_DISABLE_WARNING(warningName)
# define MSVC_DISABLE_WARNING(warningNumber)
#endif

}

#endif /* NORNIR_UTILS_HPP_ */
