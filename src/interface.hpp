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

namespace nornir{

template <typename I, typename O> class Farm;

template <typename O> class Scheduler: public AdaptiveNode{
    template <typename IN, typename OUT> friend class Farm;
private:
    ff::ff_loadbalancer* _lb;

    void setLb(ff::ff_loadbalancer* lb){
        _lb = lb;
    }
public:
    virtual ~Scheduler(){;}

    void* svc(void* task) CX11_KEYWORD(final){
        void* r = (void*) schedule();
        if(r){
           return r;
        }else{
            terminate();
            return (void*) ff::FF_EOS;
        }
    }

    virtual O* schedule() = 0;

    void send(O* task) CX11_KEYWORD(final){
        while(!ff_send_out((void*) task)){;}
    }

    void sendTo(O* task, uint id) CX11_KEYWORD(final){
        while(!_lb->ff_send_out_to(task, id)){;} //TODO: Chiedere a Massimo
    }

    void broadcast(O* task) CX11_KEYWORD(final){
        _lb->broadcast_task(task);
    }
};

template <typename I, typename O> class Worker: public AdaptiveNode{
public:
    virtual ~Worker(){;}

    void* svc(void* t) CX11_KEYWORD(final){
        return (void*) compute(reinterpret_cast<I*>(t));
    }

    virtual O* compute(I*) = 0;
};

template <typename I> class WorkerNoOut: public AdaptiveNode{
public:
    virtual ~WorkerNoOut(){;}

    void* svc(void* t) CX11_KEYWORD(final){
        compute(reinterpret_cast<I*>(t));
        return (void*) GO_ON;
    }

    virtual void compute(I*) = 0;
};

template <typename I> class Gatherer: public AdaptiveNode{
public:
    virtual ~Gatherer(){;}

    void* svc(void* t) CX11_KEYWORD(final){
        gather(reinterpret_cast<I*>(t));
        return GO_ON;
    }

    virtual void gather(I*) = 0;
};

template <typename I, typename O> class Farm{
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
    Farm(const Parameters* parameters):_scheduler(NULL), _gatherer(NULL),
                                       _params(NULL), _manager(NULL){
        _nodesCreated = false;
        _params = parameters;
        _paramsCreated = false;
    }

    Farm(const std::string& paramFileName,
         const std::string& archFileName):_scheduler(NULL), _gatherer(NULL),
                                          _params(NULL), _manager(NULL){
        _nodesCreated = false;
        _params = new Parameters(paramFileName, archFileName);
        _paramsCreated = true;
    }

    ~Farm(){
        if(_paramsCreated){
            delete _params;
        }
    }

    void addScheduler(Scheduler<I>* s){
        _farm.add_emitter(s);
    }

    void addWorker(Worker<I, O>* w){
        _workers.push_back(w);
    }

    void addGatherer(Gatherer<O>* g){
        _farm.add_collector(g);
    }

    template <class S, class W>
    void start(size_t numWorkers){
        init<S, W>(numWorkers);
        start();
    }

    template <class S, class W, class G>
    void start(size_t numWorkers){
        init<S, W>(numWorkers);
        _gatherer = new G();
        addGatherer(_gatherer);
        start();
    }

    void start(){
        _farm.add_workers(_workers);
        _scheduler->setLb(_farm.getlb());
        _manager = new ManagerFarm<>(&_farm, *_params);
        _manager->start();
    }

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
