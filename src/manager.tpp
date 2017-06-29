/*
 * manager.tpp
 *
 * Created on: 23/03/2015
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

#include "./ffincs.hpp"
#include "./parameters.hpp"
#include "./predictors.hpp"
#include "./node.hpp"
#include "./utils.hpp"

#include "external/mammut/mammut/module.hpp"
#include "external/mammut/mammut/utils.hpp"
#include "external/mammut/mammut/mammut.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <stdlib.h>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_MANAGER
#define DEBUG(x) do { cerr << "[Manager] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

class Parameters;

using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::waitForStart(){
    if(_p.qSize){
        _farm->setFixedSize(true);
        // We need to multiply for the number of workers since FastFlow
        // will divide the size for the number of workers.
        _farm->setInputQueueLength(_p.qSize * _activeWorkers.size());
        _farm->setOutputQueueLength(_p.qSize * _activeWorkers.size());
    }

    DEBUG("Init pre run");
    initNodesPreRun();

    DEBUG("Going to run");
    _farm->run_then_freeze();

    DEBUG("Init post run");
    initNodesPostRun();
    DEBUG("Farm started.");
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::askForSample(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->askForSample();
    }
}

template <typename lb_t, typename gt_t>
MonitoredSample ManagerFarm<lb_t, gt_t>::getSampleResponse(){
    MonitoredSample sample;    
    uint numActiveWorkers = _activeWorkers.size();
    for(size_t i = 0; i < numActiveWorkers; i++){
        MonitoredSample tmp;
        AdaptiveNode* w = _activeWorkers.at(i);
        w->getSampleResponse(tmp, _samples->average().latency);
        sample.loadPercentage += tmp.loadPercentage;
        sample.numTasks += tmp.numTasks;
        sample.latency += tmp.latency;
        sample.bandwidth += tmp.bandwidth;
    }
    sample.loadPercentage /= numActiveWorkers;
    sample.latency /= numActiveWorkers;
    return sample;
}

template <typename lb_t, typename gt_t>
MonitoredSample ManagerFarm<lb_t, gt_t>::getSample(){
    askForSample();
    return getSampleResponse();
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::ManagerFarm(ff_farm<lb_t, gt_t>* farm,
                                     Parameters parameters):
        Manager(parameters),
        _farm(farm),
        _emitter(dynamic_cast<AdaptiveNode*>(_farm->getEmitter())),
        _collector(dynamic_cast<AdaptiveNode*>(_farm->getCollector())),
        _activeWorkers(convertWorkers(_farm->getWorkers())){
	Manager::_pid = getpid();
    Manager::_configuration = new ConfigurationFarm(_p, _samples, _emitter,
                                                   _activeWorkers,
                                                   _collector, _farm->getgt(),
                                                   &_terminated);
    lockKnobs();
    _configuration->createAllRealCombinations();
    _selector = createSelector();
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::~ManagerFarm(){
    delete _samples;
    delete _variations;
    if(_selector){
        delete _selector;
    }
    if(Manager::_configuration){
        delete Manager::_configuration;
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initNodesPreRun() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->initPreRun(&_p, NODE_TYPE_WORKER, &_terminated);
    }
    if (_emitter) {
        _emitter->initPreRun(&_p, NODE_TYPE_EMITTER, &_terminated, _farm->getlb());
    } else {
        throw runtime_error("Emitter is needed to use the manager.");
    }
    if (_collector) {
        _collector->initPreRun(&_p, NODE_TYPE_COLLECTOR, &_terminated,
                               _farm->getgt());
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initNodesPostRun() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->initPostRun();
    }
    DEBUG("initNodesPostRun: Workers done.");
    _emitter->initPostRun();
    DEBUG("initNodesPostRun: Emitter done.");
    if (_collector) {
        _collector->initPostRun();
    }
    DEBUG("initNodesPostRun: Collector done.");
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::postConfigurationManagement(){
    const KnobVirtualCoresFarm* knobWorkers = dynamic_cast<const KnobVirtualCoresFarm*>(_configuration->getKnob(KNOB_VIRTUAL_CORES));
    std::vector<AdaptiveNode*> newWorkers = knobWorkers->getActiveWorkers();
    MonitoredSample sample;

    if(_activeWorkers.size() != newWorkers.size()){
        /**
         * Since I stopped the workers after I asked for a sample, there
         * may still be tasks that have been processed but I did not count.
         * For this reason, I get them.
         * I do not need to ask since the node put it in the Q when it
         * terminated.
         */
        DEBUG("Getting spurious..");
        sample = getSampleResponse();
        updateTasksCount(sample);
        DEBUG("Spurious got.");
    }

    _activeWorkers = newWorkers;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::terminationManagement(){
    DEBUG("Terminating...wait freezing.");
    _farm->wait_freezing();
    _farm->wait();
    DEBUG("Terminated.");
}

template <typename lb_t, typename gt_t>
ulong ManagerFarm<lb_t, gt_t>::getExecutionTime(){
    return _farm->ffTime();
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::shrinkPause(){
	KnobVirtualCoresFarm* k = ((KnobVirtualCoresFarm*) _configuration->getKnob(KNOB_VIRTUAL_CORES));
	k->prepareToFreeze();
	k->freeze();
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::stretchPause(){
	KnobVirtualCoresFarm* k = ((KnobVirtualCoresFarm*) _configuration->getKnob(KNOB_VIRTUAL_CORES));
	size_t v = k->getRealValue();
	k->prepareToRun(v);
	k->run(v);
}

}
