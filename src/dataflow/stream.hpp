/*
 * stream.hpp
 *
 * Created on: 26/03/2016
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

#ifndef NORNIR_DF_STREAM_HPP_
#define NORNIR_DF_STREAM_HPP_

#include "../../src/external/Mammut/mammut/utils.hpp"
#include "../../src/external/fastflow/ff/utils.hpp"

#include <cstdlib>

namespace nornir{
namespace dataflow{

/**
 * A generic task. The task of the skeleton (and also of the macro data flow graphs) must extend this class.
 */
class StreamElem{
public:
    virtual ~StreamElem(){};
};

/**
 * A generic wrapper. T is the type contained by the wrapper.
 */
template<typename T> class Wrapper: public StreamElem{
private:
    T x;
public:
    /**
     * Constructor of the wrapper.
     * \param x The element wrapped.
     */
    inline Wrapper(T x):x(x){;}

    /**
     * Sets the element.
     * \param t The element to set.
     */
    inline void set(T t){x=t;}

    /**
     * Returns the element wrapped.
     * \return The element wrapped.
     */
    inline T get(){return x;}
};

/**
 * A generic array wrapper. T is the type of the elements of the array.
 */
template<typename T> class ArrayWrapper:public StreamElem{
private:
    T *a;
    uint dim;
public:
    /**
     * Constructor of the wrapper.
     * \param d Size of the array.
     */
    inline ArrayWrapper(uint d):a(new T[d]),dim(d){;}

    ArrayWrapper (const ArrayWrapper& aw):a(new T[aw.dim]),dim(aw.dim){
        for(uint i=0; i<dim; i++){
            a[i]=aw.a[i];
        }
    }
    /**
     * Destructor of the wrapper.
     */
    inline ~ArrayWrapper(){delete[] a;}

    /**
     * Returns the size of the array.
     * \return The size of the array.
     */
    inline uint getSize(){return dim;}

    /**
     * Returns the array.
     * \return The array of T.
     */
    inline T* get(){return a;}
    /**
     * Returns the element at position i.
     * \return The element at position i.
     */
    inline T get(int i){return a[i];}

    /**Sets the wrapped element.**/
    inline void set(T* x){a=x;}

    /**
     * Sets the element \e x at position \e i.
     * \param i The position
     * \param x The elmenent to set.
     */
    inline void set(int i,T x){a[i]=x;}
};

template<typename T> class ArrayIndexes:public StreamElem{
private:
    ArrayWrapper<T>* aw;
    int i,j;
public:
    ArrayIndexes(ArrayWrapper<T>* aw, int i, int j):aw(aw),i(i),j(j){;}

    ArrayWrapper<T>* getArray(){return aw;}

    int geti(){return i;}

    int getj(){return j;}

};

/**
 * A generic input stream.
 */
class InputStream{
public:
    virtual ~InputStream(){;}

    /**
     * Returns the next element of the stream or NULL if no elements are
     * presents and EOS is not arrived.
     * \return The next element of the stream. NULL if no elements are presents
     * and EOS is not arrived.
     */
    virtual StreamElem* next() = 0;

    /**
     * Checks if the EndOfStream is arrived.
     * \return \e False if the EndOfStream is arrived, \e true otherwise.
     */
    virtual bool hasNext() = 0;
};

typedef struct{
    double rate;
    double duration;
}Rates;

class ClockThread: public mammut::utils::Thread{
private:
    time_t& _lastSec;
    bool& _terminated;
public:
    ClockThread(time_t& lastSec, bool& terminated);

    void run();
};

class InputStreamRate: public InputStream{
public:
private:
    std::vector<Rates> _rates;
    uint32_t _currentInterval;
    time_t _lastSec;
    time_t _startTime;
    uint64_t _processedPkts;
    time_t _currIntervalStart;
    uint32_t _currBurstSize;
    ticks _def;
    ClockThread* _clockThread;
    double _clockFrequency;
    bool _terminated;
    std::vector<StreamElem*> _objects;
    uint32_t _nextObject;
    ticks _nextBurst;
    ticks _excess;
    time_t _lastStoredRateTs;

    inline ticks ticksWait(ticks nticks) {
        ticks delta;
        ticks t0 = getticks();
        do { delta = (getticks()) - t0; } while (delta < nticks);
        return delta-nticks;
    }
protected:
    /**
     * This function must be implemented in order to load the
     * objects to be produced in the stream.
     * @return A std::vector of objects.
     **/
    virtual std::vector<StreamElem*> loadObjects() = 0;

public:
    InputStreamRate(const std::string& fileName);

    ~InputStreamRate();

    StreamElem* next();

    bool hasNext();

    void init();
};

/**
 * A generic output stream.
 */
class OutputStream{
public:
    virtual ~OutputStream(){;}

    /**
     * Puts the task into the output stream.
     * \param a The task to put.
     */
    inline virtual void put(StreamElem* a){delete a;}
};

}
}
#endif /* NORNIR_DF_STREAMS_HPP_ */
