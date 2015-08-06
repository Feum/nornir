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
#include <time.h>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_NODE
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

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
    // The percentage of time that the node spent on svc().
    double loadPercentage;

    // The number of computed tasks.
    double tasksCount;

    // The average service time (in ticks).
    double latency;

    // The bandwidth of the node.
    double bandwidthTotal;

    WorkerSample():loadPercentage(0), tasksCount(0),
                   latency(0), bandwidthTotal(0){;}

    WorkerSample& operator+=(const WorkerSample& rhs){
        loadPercentage += rhs.loadPercentage;
        tasksCount += rhs.tasksCount;
        latency += rhs.latency;
        bandwidthTotal += rhs.bandwidthTotal;
        return *this;
    }
}NodeSample;

inline WorkerSample operator+(WorkerSample lhs, const WorkerSample& rhs){
    lhs += rhs;
    return lhs;
}

/*!
 * \internal
 * \class ManagementRequestType
 * \brief Possible request types that a manager can make.
 */
typedef enum{
    // Get the current sample and reset it.
    MGMT_REQ_GET_AND_RESET_SAMPLE = 0,

    // Freezes the farm.
    MGMT_REQ_FREEZE,

    // Thaws the farm.
    MGMT_REQ_THAW
}ManagementRequestType;

/**
 * A management request.
 */
typedef struct{
    ManagementRequestType type;
    void* mark;
    size_t numWorkers;
}ManagementRequest;

/*!private
 * \class adpff_node
 * \brief This class wraps a ff_node to let it reconfigurable.
 *
 * This class wraps a ff_node to let it reconfigurable.
 */
class adpff_node: public ff_node{
private:
    friend class AdaptivityManagerFarm;

    task::TasksManager* _tasksManager;
    task::ThreadHandler* _thread;
    ManagementRequest _managementRequest;
    WorkerSample _sampleResponse;
    double _ticksPerNs;
    ticks _startTicks;

    // Queue used by the manager to notify that a request is present
    // on _managementRequest.
    ff::SWSR_Ptr_Buffer _managementQ;

    // Queue used by the node to notify that a response is present
    // on _sampleResponse.
    ff::SWSR_Ptr_Buffer _responseQ;

    double ticksToSeconds(double ticks){
        return (ticks/_ticksPerNs)/NSECS_IN_SECS;
    }

    /**
     * Initializes the node. It must be called
     * when the node is already running.
     * @param mammut A Mammut handle.
     * @param ticksPerNs The number of ticks in a nanosecond.
     */
    void init(Mammut& mammut,
              double ticksPerNs){
        size_t tid = getOSThreadId();
        if(!tid){
            throw std::runtime_error("Node init() called before "
                                     "thread creation.");
        }
        _tasksManager = mammut.getInstanceTask();
        _thread = _tasksManager->getThreadHandler(getpid(), tid);
        _ticksPerNs = ticksPerNs;
    }

    /**
     * Moves this node on a specific virtual core.
     * @param vc The virtual core where this nodes must be moved.
     */
    void move(mammut::topology::VirtualCore* vc){
        _thread->move(vc);
    }

    // Sleeps for a given amount of nanoseconds
    static inline void nsleep(long ns) {
#if defined(__linux__)
        struct timespec req={0};
        time_t sec = ns/NSECS_IN_SECS;
        req.tv_sec = sec;
        req.tv_nsec = ns - sec*NSECS_IN_SECS;
        nanosleep(&req, NULL);
#else
        throw std::runtime_error("Nanosleep not supported on this OS.");
#endif
    }

    /**
     * The result of askForSample call.
     * @param sample The statistics computed since the last
     *               time 'askForSample' has been called.
     * @param strategyPolling Strategy to apply if the queue is empty.
     * @param avgLatency Current average latency of the workers (in ns).
     * @return true if the node is running, false otherwise.
     */
    bool getSampleResponse(WorkerSample& sample,
                           StrategyPolling strategyPolling,
                           double avgLatency){
        while(_responseQ.empty()){
            if(ff_node::isfrozen()){
                return false;
            }
            switch(strategyPolling){
                case STRATEGY_POLLING_SPINNING:{
                    continue;
                }break;
                case STRATEGY_POLLING_PAUSE:{
                    PAUSE();
                }break;
                case STRATEGY_POLLING_SLEEP_SMALL:{
                    nsleep(0);
                }break;
                case STRATEGY_POLLING_SLEEP_LATENCY:{
                    nsleep(avgLatency);
                }break;
            }
        }
        _responseQ.inc();
        sample = _sampleResponse;
        return true;
    }

    /**
    * Ask the node for a sample of the statistics computed since the last
    * time this method has been called.
    * The result can be retrieved with getSampleResponse call.
    */
    void askForSample(){
        _managementRequest.type = MGMT_REQ_GET_AND_RESET_SAMPLE;
        // The value pushed in the queue will not be read,
        // it could be anything except NULL.
        _managementQ.push(&_managementRequest);
    }

    /**
     * Tells the node to freeze the farm.
     */
    void freezeAll(void* mark){
        _managementRequest.type = MGMT_REQ_FREEZE;
        _managementRequest.mark = mark;
        // The value pushed in the queue will not be read, it could be
        // anything except NULL.
        _managementQ.push(&_managementRequest);
    }

    /**
     * Thaws the farm.
     */
    void thawAll(size_t numWorkers){
        _managementRequest.type = MGMT_REQ_THAW;
        _managementRequest.numWorkers = numWorkers;
        // The value pushed in the queue will not be read, it could be
        // anything except NULL.
        _managementQ.push(&_managementRequest);
    }

    /**
     * Called on the node before starting them.
     * It must be called when the node is frozen.
     */
    void prepareToRun(){
        taskcnt = 0;
        tickstot = 0;
        _startTicks = getticks();
    }

    void storeSample(){
        int dummy;
        int* dummyPtr = &dummy;
        ticks now = getticks();
        ticks totalTicks = now - _startTicks;
        _sampleResponse.loadPercentage = ((double) (tickstot) /
                                          (double) totalTicks) * 100.0;
        _sampleResponse.tasksCount = taskcnt;
        if(taskcnt){
            _sampleResponse.latency = ((double)tickstot / (double)taskcnt) /
                                      _ticksPerNs;
        }else{
            _sampleResponse.latency = 0.0;
        }
        _sampleResponse.bandwidthTotal = (double) taskcnt /
                                         ticksToSeconds(totalTicks);

        taskcnt = 0;
        tickstot = 0;
        _startTicks = now;

        _responseQ.push(dummyPtr);
    }

    void callbackIn(void *p) CX11_KEYWORD(final){
        if(!_managementQ.empty()){
            _managementQ.inc();
            switch(_managementRequest.type){
                case MGMT_REQ_GET_AND_RESET_SAMPLE:{
                    storeSample();
                }break;
                case MGMT_REQ_FREEZE:{
                    ff_loadbalancer* lb = reinterpret_cast<ff_loadbalancer*>(p);
                    lb->broadcast_task(_managementRequest.mark);

                    /** Waits for restart request from manager. **/
                    while(_managementQ.empty()){;}
                    _managementQ.inc();
                    DEBUGB(assert(_managementRequest.type == MGMT_REQ_THAW));
                    lb->thaw(true, _managementRequest.numWorkers);
                }break;
                default:{
                    throw std::runtime_error("Unexpected mgmt request.");
                }break;
            }
        }
    }

    bool ff_send_out(void * task,
                     unsigned long retry=((unsigned long)-1),
                     unsigned long ticks=(TICKS2WAIT)) CX11_KEYWORD(final){
        callbackIn(task);
        return ff_node::ff_send_out(task, retry, ticks);
    }
public:
    /**
     * Builds an adaptive node.
     */
    adpff_node():
            _tasksManager(NULL),
            _thread(NULL),
            _managementQ(1),
            _responseQ(1){
        prepareToRun();
        _managementQ.init();
        _responseQ.init();
    }

    /**
     * Destroyes this adaptive node.
     */
    ~adpff_node(){
        if(_thread){
            _tasksManager->releaseThreadHandler(_thread);
        }
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
    virtual void notifyWorkersChange(size_t oldNumWorkers,
                                     size_t newNumWorkers){;}
};

}

#endif /* ADAPTIVE_FASTFLOW_NODE_HPP_ */
