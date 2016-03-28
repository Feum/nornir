/*
 * ewc.hpp
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

#ifndef NORNIR_DF_EWC_HPP_
#define NORNIR_DF_EWC_HPP_

#include "computable.hpp"
#include "../stream.hpp"

namespace nornir{
namespace dataflow{

/**
 * This class represents a standard pattern with an Emitter, N workers and a Collector.
 * This pattern is used in Map and Reduce skeletons.
 */
class EmitterWorkerCollector: public Computable{
private:
    uint nWorkers;
    Computable *emitter,
        *worker,
        *collector;
    bool deleteAll;
public:
    /**
     * Constructor skeleton.
     * \param emitter The skeleton used to implement the map emitter\n.
     *                Emitter can be an instance of the generic emitter
     *                (MapEmitter or ReduceEmitter) or can be redefined for the
     *                purpose simply extending the \e Computable class.
     * \param worker The skeleton used to implement the map worker.\n
     *               Worker can be an instance of the generic worker
     *               (MapWorker or ReduceWorker) or can be redefined for the
     *               purpose simply extending the \e Computable class.
     * \param collector The skeleton used to implement the map collector.\n
     *                  Collector can be an instance of the generic collector
     *                  (MapCollector or ReduceCollector) or can be redefined
     *                  for the purpose simply extending the \e Computable class.
     * \param nWorkers Number of workers.
     * \param deleteAll If true, the destructor deletes the emitter, the worker
     *                  and the collector.
     */
    inline EmitterWorkerCollector(Computable* emitter, Computable* worker,
                           Computable* collector, uint nWorkers,
                           bool deleteAll=false):
                               nWorkers(nWorkers), emitter(emitter),worker(worker),
                               collector(collector),deleteAll(deleteAll){;}

    /**
     * Destructor of the EmitterWorkerCollector.
     * If deleteAll is true, deletes emitter, worker and collector.
     */
    inline ~EmitterWorkerCollector(){
        if(deleteAll){
            delete emitter;
            delete worker;
            delete collector;
        }
    }

    /**
     * \param t Array of Task*.
     * \return Array of Task* as result. The array of Task* and its elements
     *        must be deallocated using delete[] and delete.
     *
     * This method computes the result sequentially.
    */
    inline StreamElem** compute(StreamElem** t){
        StreamElem **emitterResult,**fromWorker,
            **result=new StreamElem*[nWorkers],*toWorker[1];
        emitterResult=emitter->compute(t);

        for(uint i=0; i<nWorkers; i++){
            toWorker[0]=emitterResult[i];
            fromWorker=worker->compute(toWorker);
            result[i]=fromWorker[0];
        }
        delete[] emitterResult;
        StreamElem** fromCollector=collector->compute(result);
        delete[] result;
        return fromCollector;
    }

    /**
     * Returns the number of workers.
     * \return Number of workers.
     */
    inline uint getNWorkers(){
        return nWorkers;
    }

    /**
     * Returns a pointer to the emitter.
     * \return A pointer to the emitter.
     */
    inline Computable* getEmitter(){
        return emitter;
    }

    /**
     * Returns a pointer to the worker.
     * \return A pointer to the worker.
     */
    inline Computable* getWorker(){
        return worker;
    }

    /**
     * Returns a pointer to the collector.
     * \return A pointer to the collector.
     */
    inline Computable* getCollector(){
        return collector;
    }
};

/**
 * A generic emitter of a reduce.
 * \tparam T* is the type of the elements of the array (The input stream must produce ArrayWrapper<T*>).\n
 *         T must be a subclass of Task.
 */
template <typename T> class ReduceEmitter:public Computable{
private:
    uint chunkNum;
    bool autoDelete;
public:
    /**
     * Constructor of the emitter of the reduce.
     * \param cn Number of workers of the reduce.
     * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted.
     */
    ReduceEmitter(int cn, bool autoDelete=true);

    StreamElem** compute(StreamElem** t);
};

/**
 * A generic worker of a reduce.
 *
 * \tparam T* is the type of the elements of the array.
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 */
template <typename T, T*(*fun)(T*,T*)> class ReduceWorker:public Computable{
public:
    StreamElem** compute(StreamElem** t);
};


/**
 * A generic collector of the reduce.
 *
 * \tparam T* is the type of the elements of the array.
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 * The map returns a T*.
 */
template <typename T, T*(*fun)(T*,T*)> class ReduceCollector:public Computable{
private:
    int chunkNum;
public:
    /**
     * Constructor of the collector of the reduce.
     * \param cn Number of workers of the reduce.
     */
    ReduceCollector(int cn);

    StreamElem** compute(StreamElem** t);
};

/**
 * A generic emitter of a map.
 *
 * \tparam T* is the type of the elements of the array.(The input stream must produce ArrayWrapper<T*>).
 */
template <typename T> class MapEmitter: public Computable{
private:
    uint numWorkers;
    bool autoDelete;
public:
    /**
     * Constructor of the emitter of the map.
     * \param numWorkers Number of workers of the map.
     * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted.
     */
    MapEmitter(uint numWorkers, bool autoDelete=true);

    StreamElem** compute(StreamElem** t);
};

/**
 * Generic worker of a map.
 *
 * \tparam T* is the type of the input elements. T must be subclass of Task.
 * \tparam V* is the type of the output elements. V must be subclass of Task.
 * \tparam fun is the function to compute over the elements.
 *
 */
template <typename T, typename V, V*(*fun)(T*) > class MapWorker: public Computable{
public:
    StreamElem** compute(StreamElem** t);
};

/**
 * Generic collector of a map.
 *
 * \tparam V* is the type of the elements received from the workers.
 *
 * The \e map returns an ArrayWrapper<V*>
 */
template <typename V> class MapCollector: public Computable{
private:
    uint nWorkers;
public:
    MapCollector(uint nWorkers);

    StreamElem** compute(StreamElem** t);
};

/**
 * Creates a standard reduce.
 *
 * \tparam T is the type of the elements of the array (The input stream must produce ArrayWrapper<T*>)
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 *
 * \param nWorkers Number of workers of the reduce.
 * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted by the emitter.
 *
 * \return A pointer to a standard reduce.
 */

template<typename T, T*(*fun)(T*,T*)> EmitterWorkerCollector* createStandardReduce(int nWorkers, bool autoDelete=true){
    return new EmitterWorkerCollector(new ReduceEmitter<T>(nWorkers,autoDelete),new ReduceWorker<T,fun>,
            new ReduceCollector<T,fun>(nWorkers),nWorkers,true);
}

/**
 * Creates a standard map.
 *
 * \tparam T is the type of the input elements (The input stream must produce ArrayWrapper<T*>).
 * \tparam V is the type of the output elements (The map returns ArrayWrapper<V*>).
 * \tparam fun is the function to compute over the elements.
 *
 * \param nWorkers Number of workers of the reduce.
 * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted by the emitter.
 *
 * Input array (single element of the input stream): [x1,x2,...,xn] where typeof(x1)==typeof(x2)==...==typeof(xn)==T\n
 * Output array (single element of the output stream): [y1=fun(x1),y2=fun(x2),...,yn=fun(xn)] where typeof(y1)==typeof(y2)==...==typeof(yn)==V\n
 *
 *
 * \return A pointer to a standard map.
 */

template<typename T, typename V, V*(*fun)(T*)> EmitterWorkerCollector* createStandardMap(int nWorkers, bool autoDelete=true){
    return new EmitterWorkerCollector(new MapEmitter<T>(nWorkers,autoDelete),new MapWorker<T,V,fun>,new MapCollector<T>(nWorkers),nWorkers,true);
}

}
}

#include "ewc.tpp"

#endif /* NORNIR_DF_EWC_HPP_ */
