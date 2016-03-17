/*
 * interface.hpp
 *
 * Created on: 27/02/2016
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

#ifndef NORNIR_INTERFACE_HPP_
#define NORNIR_INTERFACE_HPP_

#include "manager.hpp"

#include <cstddef>

namespace nornir{

template <typename I, typename O> class Farm;

/**
 * @class Scheduler
 * @brief The scheduler of the farm.
 * @tparam I The type of the tasks produced by the scheduler.
 */
template <typename I> class Scheduler: public AdaptiveNode{
    template <typename IN, typename OUT> friend class Farm;
private:
    ff::ff_loadbalancer* _lb;

    void setLb(ff::ff_loadbalancer* lb){
        _lb = lb;
    }

    void* svc(void* task) CX11_KEYWORD(final){
        void* r = (void*) schedule();
        if(r){
           return r;
        }else{
            terminate();
            return (void*) ff::FF_EOS;
        }
    }
public:
    virtual ~Scheduler(){;}

    virtual I* schedule() = 0;

    /**
     * Sends a task to one of the workers.
     * @param task The task to be sent.
     */
    void send(I* task) CX11_KEYWORD(final){
        while(!ff_send_out((void*) task)){;}
    }

    /**
     * Sends a task to a specific worker.
     * @param task The task to be sent.
     * @param id The identifier of the worker. Starts from 0.
     */
    void sendTo(I* task, uint id) CX11_KEYWORD(final){
        if(id >= _lb->getnworkers()){
            throw new std::runtime_error("FATAL: Trying to send to a non running worker."
                                         "Please ensure that your application is "
                                         "notified when a rethreading occur.");
        }
        while(!_lb->ff_send_out_to(task, id)){;} //TODO: Chiedere a Massimo
    }

    /**
     * Sends a task to all the running workers.
     * @param task The task to be sent.
     */
    void broadcast(I* task) CX11_KEYWORD(final){
        _lb->broadcast_task(task);
    }
};

/**
 * @class Worker
 * @brief The worker of the farm.
 * @tparam I The type of the tasks received from the scheduler.
 * @tparam O The type of the tasks produced from the worker and sent to the
 *           gatherer. If not present, it means that the worker does not produce
 *           data.
 */
template <typename I, typename O = std::nullptr_t> class Worker: public AdaptiveNode{
private:
    void* svc(void* t) CX11_KEYWORD(final){
        return (void*) compute(reinterpret_cast<I*>(t));
    }
public:
    virtual ~Worker(){;}

    /**
     * Computes a function over a received task.
     * @param task The task received from the scheduler.
     * @return The processed task. It will be sent to the gatherer.
     */
    virtual O* compute(I* task) = 0;
};


//! @cond
/**
 * Template specialisation for the worker with no output.
 */
template <typename I>
class Worker<I, std::nullptr_t>{
private:
    void* svc(void* t) CX11_KEYWORD(final){
        compute(reinterpret_cast<I*>(t));
        return (void*) GO_ON;
    }
public:
    virtual ~Worker(){;}

    /**
     * Computes a function over a received task.
     * No output will be produced.
     * @param task The task received from the scheduler.
     */
    virtual void compute(I* task) = 0;
};
//! @endcond

/**
 * @class Gatherer
 * @brief The gatherer of the farm.
 * @tparam O The type of the tasks received from the workers.
 */
template <typename O> class Gatherer: public AdaptiveNode{
private:
    void* svc(void* t) CX11_KEYWORD(final){
        gather(reinterpret_cast<O*>(t));
        return GO_ON;
    }
public:
    virtual ~Gatherer(){;}

    //TODO: Receivefrom?

    /**
     * Computes a function over a task received from a worker.
     * @param task The task received from a worker.
     */
    virtual void gather(O* task) = 0;
};

/**
 * @class Farm
 * @brief The nornir farm.
 * @tparam I The type of the tasks produced by the scheduler and received
 *           by the workers.
 * @tparam O The type of the tasks produced by the workers and received by the
 *           gatherer
 */
template <typename I, typename O = std::nullptr_t> class Farm{
private:
    Scheduler<I>* _scheduler;
    std::vector<ff::ff_node*> _workers;
    Gatherer<O>* _gatherer;
    bool _nodesCreated;
    const Parameters* _params;
    bool _paramsCreated;
    ManagerFarm<>* _manager;
    ff::ff_farm<> _farm;

    template <class S, class W>
    void init(size_t numWorkers){
        _nodesCreated = true;
        _scheduler = new S();
        addScheduler(_scheduler);
        for(size_t i = 0; i < numWorkers; i++){
            addWorker(new W());
        }
    }
public:
    /**
     * The constructor of the farm.
     * @param The configuration parameters.
     */
    Farm(const Parameters* parameters):_scheduler(NULL), _gatherer(NULL),
                                       _params(NULL), _manager(NULL){
        _nodesCreated = false;
        _params = parameters;
        _paramsCreated = false;
    }

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    Farm(const std::string& paramFileName,
         const std::string& archFileName):_scheduler(NULL), _gatherer(NULL),
                                          _params(NULL), _manager(NULL){
        _nodesCreated = false;
        _params = new Parameters(paramFileName, archFileName);
        _paramsCreated = true;
    }

    /**
     * The destructor of the farm.
     */
    ~Farm(){
        if(_paramsCreated){
            delete _params;
        }
    }

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void addScheduler(Scheduler<I>* s){
        _farm.add_emitter(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void addWorker(Worker<I, O>* w){
        _workers.push_back(dynamic_cast<ff::ff_node*>(w));
    }

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O>* g){
        _farm.add_collector(g);
    }


    /**
     * Starts the farm.
     * @tparam S The type of the scheduler. The scheduler must have a
     *           constructor with no parameters.
     * @tparam W The type of the workers. The worker must have a constructor with
     *           no parameters.
     * @param numWorkers The maximum number of workers to be used.
     */
    template <class S, class W>
    void start(size_t numWorkers){
        init<S, W>(numWorkers);
        start();
    }

    /**
     * Starts the farm.
     * @tparam S The type of the scheduler. The scheduler must have a
     *           constructor with no parameters.
     * @tparam W The type of the workers. The worker must have a constructor with
     *           no parameters.
     * @tparam G The type of the gatherer. The gatherer must have a constructor
     *           with no parameters.
     * @param numWorkers The number of workers.
     */
    template <class S, class W, class G>
    void start(size_t numWorkers){
        init<S, W>(numWorkers);
        _gatherer = new G();
        addGatherer(_gatherer);
        start();
    }

    /**
     * Starts the farm.
     */
    void start(){
        _farm.add_workers(_workers);
        _scheduler->setLb(_farm.getlb());
        _manager = new ManagerFarm<>(&_farm, *_params);
        _manager->start();
    }

    /**
     * Waits for the farm end.
     */
    void wait(){
        if(_manager){
            _manager->join();
            delete _manager;
        }
        if(_nodesCreated){
            if(_scheduler){
                delete _scheduler;
            }
            for(size_t i = 0; i < _workers.size(); i++){
                delete _workers.at(i);
            }
            if(_gatherer){
                delete _gatherer;
            }
        }
    }
};

}

#endif /* NORNIR_INTERFACE_HPP_ */
