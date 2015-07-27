/*
 * node.hpp
 *
 * Created on: 23/03/2015
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

#ifndef ADAPTIVE_FASTFLOW_NODE_HPP_
#define ADAPTIVE_FASTFLOW_NODE_HPP_

#include <ff/node.hpp>
#include <mammut/mammut.hpp>
#include <mammut/module.hpp>
#include <mammut/utils.hpp>

#include <fstream>
#include <streambuf>
#include <string>

namespace adpff{

using namespace ff;
using namespace mammut;

/*!
 * \internal
 * \struct NodeSample
 * \brief Represents a sample of values taken from an adaptive node.
 *
 * This struct represents a sample of values taken from an adaptive node.
 */
typedef struct WorkerSample{
    double loadPercentage; ///< The percentage of time that the node spent on svc().
    double tasksCount; ///< The number of computed tasks.
    double bandwidth; ///< Bandwidth (tasks per ticks).
    double serviceTime; ///< The average service time (in ticks).
    WorkerSample():loadPercentage(0), tasksCount(0), bandwidth(0), serviceTime(0){;}

    WorkerSample& operator+=(const WorkerSample& rhs){
        loadPercentage += rhs.loadPercentage;
        tasksCount += rhs.tasksCount;
        bandwidth += rhs.bandwidth;
        serviceTime += rhs.serviceTime;
        return *this;
    }
}NodeSample;

inline WorkerSample operator+(WorkerSample lhs, const WorkerSample& rhs){
    lhs += rhs;
    return lhs;
}

/*!
 * \internal
 * \class ManagementRequest
 * \brief Possible requests that a manager can make.
 */
typedef enum{
    // Get the current sample and reset it.
    MANAGEMENT_REQUEST_GET_AND_RESET_SAMPLE = 0,

    // Produce a NULL value on output stream.
    MANAGEMENT_REQUEST_PRODUCE_NULL
}ManagementRequest;

/*!private
 * \class adpff_node
 * \brief This class wraps a ff_node to let it reconfigurable.
 *
 * This class wraps a ff_node to let it reconfigurable.
 */
class adp_ff_node: public ff_node{
private:
    template<typename lb_t, typename gt_t>
    friend class adp_ff_farm;

    friend class AdaptivityManagerFarm;

    Mammut _mammut;
    task::TasksManager* _tasksManager;
    task::ThreadHandler* _thread;
    bool _threadCreationPerformed;
    utils::Monitor _threadCreated;
    bool _threadRunning;
    utils::LockPthreadMutex _threadRunningLock;
    uint64_t _tasksCount;
    ticks _workTicks;
    ticks _startTicks;
    ManagementRequest _managementRequest;
    WorkerSample _sampleResponse;

    // Queue used by the manager to notify that a request is present
    // on _managementRequest.
    ff::SWSR_Ptr_Buffer _managementQ;

    // Queue used by the node to notify that a response is present
    // on _sampleResponse.
    ff::SWSR_Ptr_Buffer _responseQ;

    /**
     * Waits for the thread to be created.
     */
    void waitThreadCreation(){
        _threadCreated.wait();
    }

    /**
     * Returns the thread handler associated to this node.
     * If it is called before running the node, an exception
     * is thrown.
     * @return The thread handler associated to this node.
     *         It doesn't need to be released.
     */
    task::ThreadHandler* getThreadHandler() const{
        if(_thread){
            return _thread;
        }else{
            throw std::runtime_error("AdaptiveNode: Thread not initialized.");
        }
    }

    /**
     * Initializes the mammut modules needed by the node.
     * @param communicator A communicator. If NULL, the modules
     *        will be initialized locally.
     */
    void initMammutModules(Mammut& mammut){
        _mammut = mammut;
        _tasksManager = _mammut.getInstanceTask();
    }

    /**
     * Check if the node is running.
     * @return true if the node is running, false otherwise.
     */
    bool isRunning(){
        bool r;
        _threadRunningLock.lock();
        r = _threadRunning;
        _threadRunningLock.unlock();
        return r;
    }

    /**
     * Ask the node for a sample of the statistics computed since the last
     * time this method has been called.
     * The result can be retrieved with getSampleResponse call.
     */
    void askForSample(){
        _managementRequest = MANAGEMENT_REQUEST_GET_AND_RESET_SAMPLE;
        // The value pushed in the queue will not be read,
        // it could be anything except NULL.
        _managementQ.push(&_managementRequest);
    }

    /**
     * The result of askForSample call.
     * @param sample The statistics computed since the last
     *               time 'askForSample' has been called.
     * @return true if the node is running, false otherwise.
     */
     bool getSampleResponse(WorkerSample& sample){
         while(_responseQ.empty()){
             if(!isRunning()){
                 return false;
             }
         }
         _responseQ.inc();
         sample = _sampleResponse;
         return true;
     }

    /**
     * Tell the node to produce a Null task as the next task.
     */
    void produceNull(){
        _managementRequest = MANAGEMENT_REQUEST_PRODUCE_NULL;
        // The value pushed in the queue will not be read, it could be
        // anything except NULL.
        _managementQ.push(&_managementRequest);
    }

    void storeSample(){
        int dummy;
        int* dummyPtr = &dummy;
        ticks now = getticks();
        ticks totalTicks = now - _startTicks;
        _sampleResponse.loadPercentage = ((double) (_workTicks) /
                                          (double) totalTicks) * 100.0;
        _sampleResponse.tasksCount = _tasksCount;
        _sampleResponse.bandwidth = _tasksCount / (double) totalTicks;
        _sampleResponse.serviceTime = (double)_workTicks / (double)_tasksCount;

        _tasksCount = 0;
        _workTicks = 0;
        _startTicks = now;

        _responseQ.push(dummyPtr);
    }

    void resetCounters(){
        _tasksCount = 0;
        _workTicks = 0;
        _startTicks = getticks();
    }
public:
    /**
     * Builds an adaptive node.
     */
    adp_ff_node():
            _tasksManager(NULL),
            _thread(NULL),
            _threadCreationPerformed(false),
            _threadRunning(false),
            _managementRequest(MANAGEMENT_REQUEST_GET_AND_RESET_SAMPLE),
            _managementQ(1),
            _responseQ(1){
        resetCounters();
        _managementQ.init();
        _responseQ.init();
    }

    /**
     * Destroyes this adaptive node.
     */
    ~adp_ff_node(){
        if(_thread){
            _tasksManager->releaseThreadHandler(_thread);
        }
    }

    /**
     * The class that extends AdaptiveNode, must replace
     * (if present) the declaration of svc_init with
     * adp_svc_init.
     * @return 0 for success, != 0 otherwise.
     */
    virtual int adp_svc_init(){return 0;}

    /**
     * The class that extends AdaptiveNode, must replace
     * the declaration of svc with adp_svc.
     */
    virtual void* adp_svc(void* task) = 0;

    /**
     * The class that extends AdaptiveNode, must replace
     * (if present) the declaration of svc_end with
     * adp_svc_end.
     */
    virtual void adp_svc_end(){;}

    /**
     * \internal
     * Wraps the user svc_init with adaptivity logics.
     * @return 0 for success, != 0 otherwise.
     */
    int svc_init() CX11_KEYWORD(final){
        _threadRunningLock.lock();
        _threadRunning = true;
        _threadRunningLock.unlock();
        if(!_threadCreationPerformed){
            // Operations performed only the first time the thread is running.
            if(_tasksManager){
                _thread = _tasksManager->getThreadHandler();
            }else{
                throw std::runtime_error("AdaptiveNode: Tasks manager "
                                         "not initialized.");
            }
            _threadCreated.notifyAll();
            _threadCreationPerformed = true;
        }
        resetCounters();
        return adp_svc_init();
    }

    /**
     * \internal
     * Wraps the user svc with adaptivity logics.
     * @param task The input task.
     * @return The output task.
     */
    void* svc(void* task) CX11_KEYWORD(final){
        ticks start = getticks();
        void* t = adp_svc(task);
        ++_tasksCount;
        _workTicks += getticks() - start;
        assert(ff_send_out(t));

        if(!_managementQ.empty()){
            _managementQ.inc();
            switch(_managementRequest){
                case MANAGEMENT_REQUEST_GET_AND_RESET_SAMPLE:{
                    storeSample();
                }break;
                case MANAGEMENT_REQUEST_PRODUCE_NULL:{
                    return NULL;
                }
            }
        }
        return GO_ON;
    }

    /**
     * \internal
     * Wraps the user svc_end with adaptivity logics.
     */
    void svc_end() CX11_KEYWORD(final){
        storeSample();
        adp_svc_end();
        _threadRunningLock.lock();
        _threadRunning = false;
        _threadRunningLock.unlock();
    }

    /**
     * This method can be implemented by the nodes to be aware of a change
     * in the number of workers.
     * When the farm is stopped and before running it again with the new
     * number of workers, this method is called. It is called on the emitter
     * (if present), on the collector (if present) and on all the workers of
     * the new configuration.
     * In this way, if needed action may be taken to prepare for the new
     * configuration (e.g. shared state modification, etc..).
     * @param oldNumWorkers The old number of workers.
     * @param newNumWorkers The new number of workers.
     */
    virtual void notifyWorkersChange(size_t oldNumWorkers, size_t newNumWorkers){;}
};

}

#endif /* ADAPTIVE_FASTFLOW_NODE_HPP_ */
