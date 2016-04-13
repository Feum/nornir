/*
 * interface.hpp
 *
 * Created on: 27/02/2016
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

#ifndef NORNIR_INTERFACE_HPP_
#define NORNIR_INTERFACE_HPP_

#include "manager.hpp"

#include <cstddef>

namespace nornir{

template <typename I, typename O> class FarmBase;

/**
 * @class SchedulerBase
 * Common base class for scheduler (e.g. scheduler, accelerator scheduler, etc...).
 */
template <typename O> class SchedulerBase: public AdaptiveNode{
    template <typename IN, typename OUT> friend class FarmBase;
private:
    ff::ff_loadbalancer* _lb;

    void setLb(ff::ff_loadbalancer* lb){
        _lb = lb;
    }
public:
    /**
     * Sends a task to one of the workers.
     * @param task The task to be sent.
     */
    void send(O* task) CX11_KEYWORD(final){
        while(!ff_send_out((void*) task)){;}
    }

    /**
     * Sends a task to one of the workers (non blocking).
     * @param task The task to be sent.
     * @return true if the task has been sent, false otherwise.
     */
    bool sendNonBlocking(O* task) CX11_KEYWORD(final){
        return ff_send_out((void*) task);
    }

    /**
     * Sends a task to a specific worker.
     * @param task The task to be sent.
     * @param id The identifier of the worker. Starts from 0.
     */
    void sendTo(O* task, uint id) CX11_KEYWORD(final){
        if(id >= _lb->getnworkers()){
            throw new std::runtime_error("FATAL: Trying to send to a non running worker."
                                         "Please ensure that your application is "
                                         "notified when a rethreading occur.");
        }
        while(!_lb->ff_send_out_to(task, id)){;}
    }

    /**
     * Sends a task to a specific worker (non blocking).
     * @param task The task to be sent.
     * @param id The identifier of the worker. Starts from 0.
     * @return true if the task has been sent, false otherwise.
     */
    void sendToNonBlocking(O* task, uint id) CX11_KEYWORD(final){
        if(id >= _lb->getnworkers()){
            throw new std::runtime_error("FATAL: Trying to send to a non running worker."
                                         "Please ensure that your application is "
                                         "notified when a rethreading occur.");
        }
        return _lb->ff_send_out_to(task, id);
    }

    /**
     * Sends a task to all the running workers.
     * @param task The task to be sent.
     */
    void broadcast(O* task) CX11_KEYWORD(final){
        _lb->broadcast_task(task);
    }
};

/**
 * @class Scheduler
 * @brief The scheduler of the farm (accelerator case).
 * @tparam O The type of the tasks produced by the scheduler.
 */
template <typename I, typename O = std::nullptr_t> class Scheduler: public SchedulerBase<O>{
private:
    void* svc(void* task) CX11_KEYWORD(final){
        O* r = schedule((I*)(task));
        void* outTask = (void*)(r);
        if(outTask){
           return outTask;
        }else{
            AdaptiveNode::terminate();
            return (void*) ff::FF_EOS;
        }
    }
public:
    virtual ~Scheduler(){;}

    virtual O* schedule(I* task) = 0;
};


/**
 * @class Scheduler
 * @brief The scheduler of the farm.
 * @tparam O The type of the tasks produced by the scheduler.
 */
/** ATTENTION: Semantic of targs inverted!!! This is a dirty trick
 *             but we can avoid to have different class names according
 *             to the presence of an input channel.
 *             This a partial template specialisation. From the standard:
 *               The members of the class template partial specialization are
 *               unrelated to the members of the primary template. Class
 *               template partial specialization members that are used in a way
 *               that requires a definition shall be defined; the definitions of
 *               members of the primary template are never used as definitions
 *               for members of a class template partial specialization. [...]
 */
template <typename O>
class Scheduler<O, std::nullptr_t>: public SchedulerBase<O>{
private:
    void* svc(void* task) CX11_KEYWORD(final){
        O* r = schedule();
        void* outTask = (void*)(r);
        if(outTask){
           return outTask;
        }else{
            AdaptiveNode::terminate();
            return (void*) ff::FF_EOS;
        }
    }
public:
    virtual ~Scheduler(){;}

    virtual O* schedule() = 0;
};



//! @cond
/**
 * @class WorkerBase
 * Common base class for workers (e.g. worker, worker with no output, etc...).
 */
template <typename I, typename O>
class WorkerBase: public AdaptiveNode{
public:
    virtual ~WorkerBase(){;}

    /**
     * Returns the identifier of this worker.
     * @return The identifier of this worker.
     */
    uint getId() const{
        return get_my_id();
    }
};
//! @endcond

/**
 * @class Worker
 * @brief The worker of the farm.
 * @tparam I The type of the tasks received from the scheduler.
 * @tparam O The type of the tasks produced from the worker and sent to the
 *           gatherer. If not present, it means that the worker does not produce
 *           data.
 */
template <typename I, typename O = std::nullptr_t> class Worker: public WorkerBase<I, O>{
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
class Worker<I, std::nullptr_t>: public WorkerBase<I, std::nullptr_t>{
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


//! @cond
/**
 * @class GathererBase
 * @tparam I The type of data received from the workers.
 * Common base class for gatherers (e.g. gatherer, accelerator gatherer, etc...).
 */
template <typename I>
class GathererBase: public AdaptiveNode{
    //TODO: Receivefrom?
};
//! @endcond


/**
 * @class Gatherer
 * @brief The gatherer of the farm.
 * @tparam I The type of the tasks received from the workers.
 */
template <typename I, typename O = std::nullptr_t> class Gatherer: public GathererBase<I>{
private:
    void* svc(void* t) CX11_KEYWORD(final){
        O* r = gather(reinterpret_cast<I*>(t));
        return (void*)(r);
    }
public:
    virtual ~Gatherer(){;}

    /**
     * Computes a function over a task received from a worker.
     * @param task The task received from a worker.
     */
    virtual O* gather(I* task) = 0;
};

/**
 * Template specialisation for the worker with no output.
 */
template <typename I>
class Gatherer<I, std::nullptr_t>: public GathererBase<I>{
private:
    void* svc(void* t) CX11_KEYWORD(final){
        gather(reinterpret_cast<I*>(t));
        return GO_ON;
    }
public:
    virtual ~Gatherer(){;}

    /**
     * Computes a function over a task received from a worker.
     * @param task The task received from a worker.
     */
    virtual void gather(I* task) = 0;
};


/**
 * @class FarmBase
 * @brief The nornir farm.
 * @tparam I The type of the tasks produced by the scheduler and received
 *           by the workers.
 * @tparam O The type of the tasks produced by the workers and received by the
 *           gatherer
 */
template <typename I, typename O> class FarmBase{
protected:
    ff::ff_farm<>* _farm;
    std::vector<ff::ff_node*> _workers;
private:
    bool _nodesCreated;
    const Parameters* _params;
    bool _paramsCreated;
    ManagerFarm<>* _manager;

    template <class S, class W>
    void init(size_t numWorkers){
        _nodesCreated = true;
        setScheduler(new S());
        for(size_t i = 0; i < numWorkers; i++){
            setWorker(new W());
        }
    }

protected:
    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void setScheduler(SchedulerBase<I>* s){
        if(!s){
            throw std::runtime_error("setScheduler: Scheduler must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _farm->add_emitter(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void setWorker(WorkerBase<I, O>* w){
        if(!w){
            throw std::runtime_error("setWorker: Worker must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _workers.push_back(dynamic_cast<ff::ff_node*>(w));
    }

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void setGatherer(GathererBase<O>* g){
        if(!g){
            throw std::runtime_error("setGatherer: Gatherer must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _farm->add_collector(g);
    }

    virtual void createFarm(){
        _farm = new ff::ff_farm<>();
    }

    virtual void preStart(){
        ;
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    FarmBase(const Parameters* parameters):_farm(NULL),
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
    FarmBase(const std::string& paramFileName,
         const std::string& archFileName):_farm(NULL),
                                          _params(NULL), _manager(NULL){
        _nodesCreated = false;
        _params = new Parameters(paramFileName, archFileName);
        _paramsCreated = true;
    }

    /**
     * The destructor of the farm.
     */
    virtual ~FarmBase(){
        if(_paramsCreated){
            delete _params;
        }

        if(_farm){
            delete _farm;
        }
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
        setGatherer(new G());
        start();
    }

    /**
     * Starts the farm.
     */
    void start(){
        if(!_farm){
            createFarm();
        }
        _farm->add_workers(_workers);
        if(_farm->getEmitter()){
            (dynamic_cast<SchedulerBase<I>*>(_farm->getEmitter()))->setLb(_farm->getlb());
        }
        preStart();
        _manager = new ManagerFarm<>(_farm, *_params);
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
            if(_farm->getEmitter()){
                delete _farm->getEmitter();
            }
            for(size_t i = 0; i < _workers.size(); i++){
                delete _workers.at(i);
            }
            if(_farm->getCollector()){
                delete _farm->getCollector();
            }
        }
    }

    /**
     * Returns the number of workers currently present in the farm.
     */
    size_t getCurrentNumWorkers() const{
        return _farm->getlb()->getnworkers();
    }

    void stats(std::ostream& o) const{
        _farm->ffStats(o);
    }
};

/**
 * @class Farm
 * @brief The nornir farm.
 * @tparam I The type of the tasks produced by the scheduler and received
 *           by the workers.
 * @tparam O The type of the tasks produced by the workers and received by the
 *           gatherer
 */
template <typename I, typename O = std::nullptr_t> class Farm: public FarmBase<I, O>{
public:
    /**
     * The constructor of the farm.
     * @param parameter The configuration parameters.
     */
    Farm(const Parameters* parameters):FarmBase<I,O>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    Farm(const std::string& paramFileName,
         const std::string& archFileName):FarmBase<I,O>(paramFileName, archFileName){;}

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void addScheduler(Scheduler<I>* s){
        FarmBase<I,O>::setScheduler(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void addWorker(Worker<I, O>* w){
        FarmBase<I,O>::setWorker(w);
    }

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O>* g){
        FarmBase<I,O>::setGatherer(g);
    }
};


template <typename I, typename O> class SchedulerDummy: public Scheduler<I, O>{
public:
    O* schedule(I* task){
        return (O*) task;
    }
};

template <typename I> class GathererDummy: public Gatherer<I, I>{
public:
    I* gather(I* task){
        return task;
    }
};

/**
 * @class FarmAcceleratorBase
 * @brief The base class for farm accelerators.
 * @tparam S The type of data sent from the application to the scheduler.
 * @tparam I The type of data sent by the scheduler to the workers.
 * @tparam O The type of data sent by the workers to the gatherer.
 * @tparam G The type of data sent by the gatherer to the application.
 */
template <typename S, typename I, typename O, typename G>
class FarmAcceleratorBase: public FarmBase<I, O>{
private:
    SchedulerDummy<S, I>* _schedulerDummy;
protected:
    void createFarm(){
        FarmBase<I, O>::_farm = new ff::ff_farm<>(true);
    }

    void preStart(){
        // Adds default gatherer and scheduler.
        if(!FarmBase<I, O>::_farm->getEmitter()){
            _schedulerDummy = new SchedulerDummy<S, I>();
            FarmBase<I, O>::setScheduler(_schedulerDummy);
        }
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    FarmAcceleratorBase(const Parameters* parameters):
            FarmBase<I,O>::FarmBase(parameters),
            _schedulerDummy(NULL){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    FarmAcceleratorBase(const std::string& paramFileName,
                    const std::string& archFileName):
             FarmBase<I,O>::FarmBase(paramFileName, archFileName),
             _schedulerDummy(NULL){;}

    /**
     * Denstructor of the accelerator.
     */
    ~FarmAcceleratorBase(){
        if(_schedulerDummy){
            delete _schedulerDummy;
        }
    }

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void addScheduler(Scheduler<S, I>* s){
        FarmBase<I,O>::setScheduler(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void addWorker(Worker<I, O>* w){
        FarmBase<I,O>::setWorker(w);
    }

    /**
     * Offloads a task to the accelerator.
     * If the task is NULL, the accelerator will be shutdown. You need to
     * wait for its termination with the wait() call.
     * This call is blocking and only returns when the task has been sent.
     * @param task The task to be offloaded.
     */
    void offload(S* task){
        void* realTask;
        if(!task){
            realTask = EOS;
        }else{
            realTask = (void*) task;
        }
        while(!FarmBase<I,O>::_farm->offload(realTask));
        if(realTask == EOS){
            ((Scheduler<S, I>*) FarmBase<I, O>::_farm->getEmitter())->terminate();
        }
    }

    /**
     * Offloads a task to the accelerator.
     * If the task is NULL, the accelerator will be shutdown. You need to
     * wait for its termination with the wait() call.
     * This call is non blocking.
     * @param task The task to be offloaded.
     * @return True if the task has been sent, false otherwise.
     */
    bool offloadNonBlocking(S* task){
        void* realTask;
        if(!task){
            realTask = EOS;
        }else{
            realTask = (void*) task;
        }
        bool res = FarmBase<I,O>::_farm->offload(realTask, 1);
        if(res && realTask == EOS){
            ((Scheduler<S, I>*) FarmBase<I, O>::_farm->getEmitter())->terminate();
        }
        return res;
    }

    /**
     * Shutdown the farm. You need to wait for its termination with the
     * wait() call.
     * This call blocks until it is possible to send a shutdown request.
     */
    void shutdown(){
        offload(NULL);
    }

    /**
     * Shutdown the farm. You need to wait for its termination with the
     * wait() call.
     * @return True if the shutdown request has been sent, false otherwise.
     */
    bool shutdownNonBlocking(){
        return offloadNonBlocking(NULL);
    }
};

/**
 * @class FarmAccelerator
 * @brief Accelerator with collector and output channel.
 *
 * @tparam S The type of data sent from the application to the scheduler.
 * @tparam I The type of data sent by the scheduler to the workers.
 * @tparam O The type of data sent by the workers to the gatherer.
 * @tparam G The type of data sent by the gatherer to the application.
 *
 * - 1 tparams specified: Accelerator with no collector. Same type between
 *     application and scheduler and between scheduler and workers.
 * - 2 tparams specified: Accelerator with no collector. One type between
 *     application and scheduler and the other between scheduler and workers.
 *     You MUST explicitly provide a Scheduler.
 * - 3 tparams specified: Accelerator with collector. One type between application
 *     and scheduler, one type between scheduler and workers, and one type between
 *     workers and collector.
 * - 4 tparams specified: Accelerator with collector. One type for each channel(s).
 *
 */
template <typename S, typename I = S, typename O = std::nullptr_t, typename G = std::nullptr_t>
class FarmAccelerator: public FarmAcceleratorBase<S, I, O, G>{
private:
    GathererDummy<O>* _gathererDummy;

    bool isManagementTask(void* task){
        return task == EOSW   ||
               task == GO_OUT ||
               task == BLK    ||
               task == NBLK;
    }
protected:
    void preStart(){
        FarmAcceleratorBase<S, I, O, G>::preStart();
        if(!FarmBase<I, O>::_farm->getCollector()){
            _gathererDummy = new GathererDummy<O>();
            FarmBase<I, O>::setGatherer(_gathererDummy);
        }
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, O, G>(parameters),
            _gathererDummy(NULL){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    FarmAccelerator(const std::string& paramFileName,
                    const std::string& archFileName):
            FarmAcceleratorBase<S, I, O, G>(paramFileName, archFileName),
            _gathererDummy(NULL){;}

    /**
     * Denstructor of the accelerator.
     */
    ~FarmAccelerator(){
        if(_gathererDummy){
            delete _gathererDummy;
        }
    }


    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O, G>* g){
        FarmBase<I, O>::setGatherer(g);
    }

    /**
     * Gets a result from the accelerator.
     * This call is blocking and only returns when the task has been received.
     * @return The retrieved result. If NULL, then this is the last produced
     *         result.
     */
    G* getResult(){
        G* r = NULL;
        G** pr = &r;
        while((!FarmBase<I,O>::_farm->load_result((void**)pr) && *pr != EOS) ||
              isManagementTask((void*)*pr)){
            ;
        }

        if(*pr == EOS){
            *pr = NULL;
        }
        return *pr;
    }

    /**
     * Gets a result from the accelerator.
     * This call is non blocking.
     * @param validResult Will be set to true if a result has been retrieved,
     *        to false otherwise.
     * @return The retrieved result. If NULL, then this is the last produced
     *         result.
     */
    G* getResultNonBlocking(bool& validResult){
        G* r = NULL;
        G** pr = &r;
        if(!FarmBase<I,O>::_farm->load_result_nb((void**)pr) ||
           isManagementTask((void*)*pr)){
            validResult = false;
            return NULL;
        }else{
            validResult = true;
            if(*pr == EOS){
                *pr = NULL;
            }
            return *pr;
        }
    }

};

/**
 * @brief Accelerator with no collector.
 */
template <typename S, typename I>
class FarmAccelerator<S, I, std::nullptr_t, std::nullptr_t>: public FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>{
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    FarmAccelerator(const std::string& paramFileName,
                    const std::string& archFileName):
            FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>(paramFileName, archFileName){;}
};

/**
 * @brief Accelerator with collector but no output channel.
 */
template <typename S, typename I, typename O>
class FarmAccelerator<S, I, O, std::nullptr_t>: public FarmAcceleratorBase<S, I, O, std::nullptr_t>{
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, O, std::nullptr_t>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     * @param archFileName The filename of the XML file containing the
     *        architectural parameters.
     */
    FarmAccelerator(const std::string& paramFileName,
                    const std::string& archFileName):
            FarmAcceleratorBase<S, I, O, std::nullptr_t>(paramFileName, archFileName){;}

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O>* g){
        FarmBase<I, O>::setGatherer(g);
    }
};


}

#endif /* NORNIR_INTERFACE_HPP_ */
