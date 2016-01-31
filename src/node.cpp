/*
 * node.cpp
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

#include "./node.hpp"
#include "./utils.hpp"

#include <mammut/module.hpp>
#include <mammut/utils.hpp>

#include <fstream>
#include <streambuf>
#include <string>
#include <time.h>

#include <ff/lb.hpp>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_NODE
#define DEBUG(x) do { cerr << "[Node] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace adpff{

using namespace std;

using namespace ff;
using namespace mammut;
using namespace mammut::task;
using namespace mammut::topology;


// Sleeps for a given amount of nanoseconds
static inline void nSleep(long ns) {
#if defined(__linux__)
    struct timespec req={0};
    time_t sec = ns/NSECS_IN_SECS;
    req.tv_sec = sec;
    req.tv_nsec = ns - sec*NSECS_IN_SECS;
    nanosleep(&req, NULL);
#else
    throw runtime_error("Nanosleep not supported on this OS.");
#endif
}

double AdaptiveNode::ticksToSeconds(double ticks) const{
    return (ticks/_ticksPerNs)/NSECS_IN_SECS;
}

void AdaptiveNode::initPreRun(Mammut& mammut, double ticksPerNs,
                              NodeType nodeType, volatile bool* terminated,
                              ff::ff_thread* ffThread){
    _tasksManager = mammut.getInstanceTask();
    if(!_tasksManager){
        throw runtime_error("Node init(): impossible to "
                            "get the tasks manager.");
    }
    _ticksPerNs = ticksPerNs;
    _nodeType = nodeType;
    _terminated = terminated;
    _ffThread = ffThread;
}

void AdaptiveNode::initPostRun(){
    DEBUG("Waiting for start.");
    while(!_started){;}
    DEBUG("Started.");
    size_t tid;
    if(_ffThread){
        tid = _ffThread->getOSThreadId();
    }else{
        tid = getOSThreadId();
    }
    assert(tid);

    _thread = _tasksManager->getThreadHandler(getpid(), tid);
}

void AdaptiveNode::clean(){
    if(_thread){
        _tasksManager->releaseThreadHandler(_thread);
        _thread = NULL;
    }
}

void AdaptiveNode::move(VirtualCore* vc){
    if(_thread){
        _thread->move(vc);
    }
}

void AdaptiveNode::getSampleResponse(WorkerSample& sample,
                       StrategyPolling strategyPolling,
                       double avgLatency){
    while(_responseQ.empty()){
        if(*_terminated){return;}
        switch(strategyPolling){
            case STRATEGY_POLLING_SPINNING:{
                continue;
            }break;
            case STRATEGY_POLLING_PAUSE:{
                PAUSE();
            }break;
            case STRATEGY_POLLING_SLEEP_SMALL:{
                nSleep(0);
            }break;
            case STRATEGY_POLLING_SLEEP_LATENCY:{
                nSleep(avgLatency);
            }break;
        }
    }
    _responseQ.inc();
    sample = _sampleResponse;
}

void AdaptiveNode::resetSample(){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_RESET_SAMPLE];
    request->type = MGMT_REQ_RESET_SAMPLE;
    while(!_managementQ.push(request));
    DEBUG("RESETSAMPLE");
  }


void AdaptiveNode::askForSample(){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_GET_AND_RESET_SAMPLE];
    request->type = MGMT_REQ_GET_AND_RESET_SAMPLE;
    // The value pushed in the queue will not be read,
    // it could be anything except NULL.
    while(!_managementQ.push(request));
    DEBUG("ASKFORSAMPLE");
}

void AdaptiveNode::setQBlocking(){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_SWITCH_BLOCKING];
    request->type = MGMT_REQ_SWITCH_BLOCKING;
    request->mark = (void*) FF_BLK;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    while(!_managementQ.push(request));
    DEBUG("Queues switched to blocking.");
}

void AdaptiveNode::setQNonblocking(){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_SWITCH_BLOCKING];
    request->type = MGMT_REQ_SWITCH_BLOCKING;
    request->mark = (void*) FF_NBLK;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    while(!_managementQ.push(request));
    DEBUG("Queues switched to nonblocking.");
}

void AdaptiveNode::freezeAll(void* mark){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_FREEZE];
    request->type = MGMT_REQ_FREEZE;
    request->mark = mark;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    DEBUG("FREEZEALLBEF");
    while(!_managementQ.push(request));
    DEBUG("FREEZEALLAFT");
    if(_nodeType == NODE_TYPE_EMITTER){
        DEBUG("Pushed freeze req");
    }
}

void AdaptiveNode::thawAll(size_t numWorkers){
    ManagementRequest* request = &_managementRequests[MGMT_REQ_THAW];
    request->type = MGMT_REQ_THAW;
    request->numWorkers = numWorkers;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    DEBUG("THAWALLBEF");
    while(!_managementQ.push(request));
    DEBUG("THAWALLAFT");
}

void AdaptiveNode::prepareToFreeze(){
    ;
}

void AdaptiveNode::prepareToRun(){
    reset();
}

void AdaptiveNode::reset(){
    tickstot = 0;
    taskcnt = 0;
    _tasksCount = 0;
    _ticksWork = 0;
    _startTicks = getticks();
}

void AdaptiveNode::storeSample(){
    DEBUG("Storing sample");
    int dummy;
    int* dummyPtr = &dummy;
    ticks totalTicks = getticks() - _startTicks; // Including idle periods

    _ticksWork = tickstot;
    _tasksCount = taskcnt;

    _sampleResponse.loadPercentage = ((double) (_ticksWork) / (double) totalTicks)
                                     * 100.0;
    _sampleResponse.tasksCount = _tasksCount;
    if(_tasksCount){
        _sampleResponse.latency = ((double)_ticksWork / (double)_tasksCount) /
                                  _ticksPerNs;
    }else{
        _sampleResponse.latency = 0.0;
    }
    _sampleResponse.bandwidthTotal = (double) _tasksCount / ticksToSeconds(totalTicks);

    reset();

    assert(_responseQ.push(dummyPtr));
}

void AdaptiveNode::callbackIn(void *p) CX11_KEYWORD(final){
    _started = true;

    ManagementRequest* request;

    while(_managementQ.pop((void**) &request)){
        switch(request->type){
            case MGMT_REQ_GET_AND_RESET_SAMPLE:{
                DEBUG("Get and reset received");
                storeSample();
            }break;
            case MGMT_REQ_RESET_SAMPLE:{
                DEBUG("Reset received");
                reset();
            }break;
            case MGMT_REQ_FREEZE:{
                /**
                 * When the emitter returns the EOS, a broadcast will be executed.
                 * The broadcast will call the callbacks, so we could pop
                 * a freeze request while the farm is already terminated.
                 * For this reason we must check the flag.
                 *
                 * TODO: This should be useless on the new ff support since the
                 * callback should now be called only on real task and not
                 * when marks are received (for the moment we still do the
                 * check because this is not always true).
                 */
                if(!*_terminated){
                    assert(_nodeType == NODE_TYPE_EMITTER);
                    DEBUG("Freeze request received");
                    ff_loadbalancer* lb = reinterpret_cast<ff_loadbalancer*>(p);
                    lb->broadcast_task(request->mark);
                    DEBUG("Broadcasted");
                    svc_end();

                    /** Waits for restart request from manager. **/
                    while(!_managementQ.pop((void**) &request));
                    DEBUGB(assert(request->type == MGMT_REQ_THAW));
                    svc_init();
                    lb->thawWorkers(true, request->numWorkers);
                }
            }break;
            case MGMT_REQ_SWITCH_BLOCKING:{
                assert(_nodeType == NODE_TYPE_EMITTER);
                DEBUG("Block/Nonblock request received");
                ff_loadbalancer* lb = reinterpret_cast<ff_loadbalancer*>(p);
                lb->broadcast_task(request->mark);
                DEBUG("Broadcasted");
            }break;
            default:{
                throw runtime_error("Unexpected mgmt request.");
            }break;
        }
    }
}

void AdaptiveNode::callbackOut(void *p) CX11_KEYWORD(final){
    callbackIn(p);
}

void AdaptiveNode::eosnotify(ssize_t id) CX11_KEYWORD(final){
    if(_nodeType == NODE_TYPE_WORKER){
        storeSample();
    }
}



AdaptiveNode::AdaptiveNode():
        _started(false),
        _terminated(NULL),
        _tasksManager(NULL),
        _thread(NULL),
        _ticksWork(0),
        _tasksCount(0),
        // Some messages are without an answer (e.g. SWITCH_BLOCKING or
        // RESET_SAMPLE). For this reason, we could enqueue more request before
        // the node reads any of them. For example, we could enqueue a
        // SWITCH_BLOCKING, a RESET_SAMPLE and a GET_AND_RESET_SAMPLE.
        // For this reason, the size of the managementQ is greater than 1.
        _managementQ(4),
        _responseQ(2){
    _managementQ.init();
    _responseQ.init();
    prepareToRun();
}

AdaptiveNode::~AdaptiveNode(){
    clean();
}

void AdaptiveNode::terminate(){
    *_terminated = true;
}

void AdaptiveNode::notifyWorkersChange(size_t oldNumWorkers,
                                       size_t newNumWorkers){;}
}

