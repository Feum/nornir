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

void AdaptiveNode::initPreRun(Mammut& mammut, double ticksPerNs, NodeType nodeType){
    _tasksManager = mammut.getInstanceTask();
    if(!_tasksManager){
        throw runtime_error("Node init(): impossible to "
                            "get the tasks manager.");
    }
    _ticksPerNs = ticksPerNs;
    _nodeType = nodeType;
}

void AdaptiveNode::initPostRun(){
    while(!_started){;}
    size_t tid = getOSThreadId();
    assert(tid != 0);
    _thread = _tasksManager->getThreadHandler(getpid(), tid);
    if(_nodeType == NODE_TYPE_EMITTER){
        DEBUG("EMITTERTID: " << tid);
    }
}

void AdaptiveNode::clean(){
    if(_thread){
        _tasksManager->releaseThreadHandler(_thread);
        _thread = NULL;
    }
}

void AdaptiveNode::move(VirtualCore* vc){
    _thread->move(vc);
}

void AdaptiveNode::getSampleResponse(WorkerSample& sample,
                       StrategyPolling strategyPolling,
                       double avgLatency){
    while(_responseQ.empty()){
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

void AdaptiveNode::askForSample(){
    _managementRequest.type = MGMT_REQ_GET_AND_RESET_SAMPLE;
    // The value pushed in the queue will not be read,
    // it could be anything except NULL.
    while(!_managementQ.push(&_managementRequest));
    DEBUG("ASKFORSAMPLE");
}

void AdaptiveNode::freezeAll(void* mark){
    _managementRequest.type = MGMT_REQ_FREEZE;
    _managementRequest.mark = mark;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    DEBUG("FREEZEALLBEF");
    while(!_managementQ.push(&_managementRequest));
    DEBUG("FREEZEALLAFT");
    if(_nodeType == NODE_TYPE_EMITTER){
        DEBUG("Pushed freeze req");
    }
}

void AdaptiveNode::thawAll(size_t numWorkers){
    _managementRequest.type = MGMT_REQ_THAW;
    _managementRequest.numWorkers = numWorkers;
    // The value pushed in the queue will not be read, it could be
    // anything except NULL.
    DEBUG("THAWALLBEF");
    while(!_managementQ.push(&_managementRequest));
    DEBUG("THAWALLAFT");
}

void AdaptiveNode::prepareToFreeze(){
    _goingToFreeze = true;
}

void AdaptiveNode::prepareToRun(){
    while(!_responseQ.empty()){
        DEBUG("Clearing response Q.");
        _responseQ.inc();
    }
    _startTicks = getticks();
    _goingToFreeze = false;
    _tasksCount = 0;
    _ticksTot = 0;
}

bool AdaptiveNode::isTerminated() const{
    return _terminated;
}

void AdaptiveNode::storeSample(){
    int dummy;
    int* dummyPtr = &dummy;
    ticks now = getticks();
    ticks totalTicks = now - _startTicks;
    _sampleResponse.loadPercentage = ((double) (_ticksTot) /
                                      (double) totalTicks) * 100.0;
    _sampleResponse.tasksCount = _tasksCount;
    if(_tasksCount){
        _sampleResponse.latency = ((double)_ticksTot / (double)_tasksCount) /
                                  _ticksPerNs;
    }else{
        _sampleResponse.latency = 0.0;
    }
    _sampleResponse.bandwidthTotal = (double) _tasksCount /
                                     ticksToSeconds(totalTicks);

    _tasksCount = 0;
    _ticksTot = 0;
    _startTicks = now;

    assert(_responseQ.push(dummyPtr));
}

void AdaptiveNode::callbackIn(void *p) CX11_KEYWORD(final){
    _started = true;

    if(!_managementQ.empty()){
        _managementQ.inc();
        if(_nodeType == NODE_TYPE_EMITTER){
            DEBUG("Popped something");
        }
        switch(_managementRequest.type){
            case MGMT_REQ_GET_AND_RESET_SAMPLE:{
                DEBUG("Get and reset received");
                storeSample();
            }break;
            case MGMT_REQ_FREEZE:{
                DEBUG("Freeze request received");
                ff_loadbalancer* lb = reinterpret_cast<ff_loadbalancer*>(p);
                lb->broadcast_task(_managementRequest.mark);
                DEBUG("Broadcasted");

                /** Waits for restart request from manager. **/
                while(_managementQ.empty());
                _managementQ.inc();
                DEBUGB(assert(_managementRequest.type == MGMT_REQ_THAW));
                lb->thawWorkers(true, _managementRequest.numWorkers);
            }break;
            default:{
                throw runtime_error("Unexpected mgmt request.");
            }break;
        }
    }
}

void AdaptiveNode::eosnotify(ssize_t) CX11_KEYWORD(final){
    if(!_goingToFreeze){
        _terminated = true;
        if(_nodeType == NODE_TYPE_WORKER){
            storeSample();
        }
    }
}

AdaptiveNode::AdaptiveNode():
        _started(false),
        _terminated(false),
        _goingToFreeze(false),
        _tasksManager(NULL),
        _thread(NULL),
        _ticksTot(0),
        _tasksCount(0),
        _managementQ(1),
        _responseQ(2){
    _managementQ.init();
    _responseQ.init();
    prepareToRun();
}

AdaptiveNode::~AdaptiveNode(){
    clean();
}

void AdaptiveNode::notifyWorkersChange(size_t oldNumWorkers,
                                       size_t newNumWorkers){;}
}

