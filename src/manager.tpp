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
MonitoredSample ManagerFarm<lb_t, gt_t>::getSample(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->askForSample();
    }
    MonitoredSample sample;
    AdaptiveNode* w;
    uint numActiveWorkers = _activeWorkers.size();

    for(size_t i = 0; i < numActiveWorkers; i++){
        MonitoredSample tmp;
        w = _activeWorkers.at(i);
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
    DEBUGB(samplesFile.open("samples.csv"));
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
    DEBUGB(samplesFile.close());
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
        sample = getSample();
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

typedef struct{
    double numCores;
    double frequency;
}SimulationKey;

inline std::ostream& operator<<(std::ostream& os, const SimulationKey& obj){
    os << "[";
    os << obj.numCores << ", " << obj.frequency;
    os << "]";
    return os;
}

struct SimulationKeyCompare{
   bool operator()(const SimulationKey& lhs, const SimulationKey& rhs) const{
       if(lhs.numCores != rhs.numCores){
           return lhs.numCores < rhs.numCores;
       }else{
           return lhs.frequency < rhs.frequency;
       }
   }
};

typedef struct{
    double completionTime;
    double wattsCpu;
    double wattsCores;
    double wattsDram;
}SimulationData;

template <typename lb_t, typename gt_t>
SimulationResult ManagerFarm<lb_t, gt_t>::simulate(std::vector<std::string>& configurationData, volatile bool* terminate, size_t maxIterations){
    // Deprecated since we do not have anymore primary/secondary values.
#if 0
    vector<string>& lines = configurationData;
    map<SimulationKey, SimulationData, SimulationKeyCompare> table;
    KnobsValues lastConfigurationValues = _configuration->getRealValues();
    // Starts from 1 to skip the header line.
    for(size_t i = 0; i < lines.size(); i++){
        SimulationKey key;
        SimulationData data;
        vector<string> fields = split(lines.at(i), '\t');
        key.numCores = atof(fields[0].c_str());
        key.frequency = atof(fields[1].c_str());

        data.completionTime = atof(fields[2].c_str());
        data.wattsCpu = atof(fields[3].c_str());
        if(fields.size() > 4){
            data.wattsCores = atof(fields[4].c_str());
        }
        if(fields.size() > 5){
            data.wattsDram = atof(fields[5].c_str());
        }

        table.insert(std::pair<SimulationKey, SimulationData>(key, data));
    }


    if(_p.requirements.executionTime != NORNIR_REQUIREMENT_UNDEF){
        _remainingTasks = _p.requirements.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 + _p.requirements.executionTime;
    }

    waitForStart();

    if(_counter){
        _counter->reset();
    }
    _lastStoredSampleMs = getMillisecondsTime();
    if(_p.observer){
        _p.observer->_startMonitoringMs = _lastStoredSampleMs;
    }

    /* Force the first calibration point. **/
    KnobsValues kv = decide();
    act(kv);

    _samples->reset();
    _variations->reset();

    SimulationResult res;
    KnobsValues values = _configuration->getRealValues();
    SimulationKey key;
    SimulationData data;
    key.numCores = values[KNOB_VIRTUAL_CORES];
    key.frequency = values[KNOB_FREQUENCY];
    if(table.find(key) != table.end()){
        data = table[key];
    }else{
        throw std::runtime_error("Impossible to find value for key.");
    }

    MonitoredSample sample;
    sample.watts = data.wattsCores;
    sample.bandwidth = 1.0 / data.completionTime;
    sample.loadPercentage = 100.0; //TODO Not yet implemented.
    sample.latency = 0; //TODO Not yet implemented.
    size_t steps = 1;

    if(_p.synchronousWorkers){
        sample.bandwidth /= key.numCores;
    }

    res.currentBandwidth = sample.bandwidth;
    res.currentPower = sample.watts;

    _samples->add(sample);
    _variations->add(getPrimaryValue(_samples->coefficientVariation()));


    ThreadHandler* thisThread = _task->getProcessHandler(getpid())->getThreadHandler(gettid());
    thisThread->move(MANAGER_VIRTUAL_CORE);

    double microsecsSleep = 0;
    double startSample = getMillisecondsTime();
    double overheadMs = 0;

    uint samplingInterval;
    uint steadySamples = 0;

    while((!_configuration->equal(lastConfigurationValues) || _selector->isCalibrating()) && (!maxIterations || steps <= maxIterations)){
        ++steps;
        lastConfigurationValues = _configuration->getRealValues();
        overheadMs = getMillisecondsTime() - startSample;
        if(_selector->isCalibrating()){
            samplingInterval = _p.samplingIntervalCalibration;
            steadySamples = 0;
        }else if(steadySamples < _p.steadyThreshold){
            samplingInterval = _p.samplingIntervalCalibration;
            steadySamples++;
        }else{
            samplingInterval = _p.samplingIntervalSteady;
        }
        microsecsSleep = ((double)samplingInterval - overheadMs)*
                          (double)MAMMUT_MICROSECS_IN_MILLISEC;
        if(microsecsSleep < 0){
            microsecsSleep = 0;
        }else{
            usleep(microsecsSleep);
        }

        startSample = getMillisecondsTime();
        DEBUG("Storing new sample.");
        observe();
        _samples->reset();
        _variations->reset();


        MonitoredSample sample;
        sample.watts = data.wattsCores;
        sample.bandwidth = 1.0 / data.completionTime;
        sample.loadPercentage = 100.0; //TODO Not yet implemented.
        sample.latency = 0; //TODO Not yet implemented.

        if(_p.synchronousWorkers){
            sample.bandwidth /= key.numCores;
        }

        res.currentBandwidth = sample.bandwidth;
        res.currentPower = sample.watts;

        _samples->add(sample);
        _variations->add(getPrimaryValue(_samples->coefficientVariation()));

        DEBUG("New sample stored.");

        updateRequiredBandwidth();

        logObservation();

        if(!persist()){
            DEBUG("Asking selector.");
            KnobsValues kv = decide();
            act(kv);

            values = _configuration->getRealValues();
            key.numCores = values[KNOB_VIRTUAL_CORES];
            key.frequency = values[KNOB_FREQUENCY];
            if(table.find(key) != table.end()){
                data = table[key];
            }else{
                throw std::runtime_error("Impossible to find value for key.");
            }

            _configuration->trigger();
            startSample = getMillisecondsTime();
        }
    }

    *terminate = true;
    terminationManagement();

    // We do -1 because the last step was the optimal configuration.
    res.numSteps = steps - 1;
    res.foundConfiguration = _configuration->getRealValues();

    switch(_p.strategySelection){
        // For SelectorPredictive selectors we compute the MAPE.
        case STRATEGY_SELECTION_LEARNING:
        case STRATEGY_SELECTION_MISHRA:
        case STRATEGY_SELECTION_FULLSEARCH:
        case STRATEGY_SELECTION_ANALYTICAL:{
            double primaryPrediction, secondaryPrediction;
            double primaryValue, secondaryValue;
            size_t keys = 0;
            double mapeBandwidth = 0, mapePower = 0;
            SelectorPredictive* sel = (SelectorPredictive*) _selector;
            std::vector<KnobsValues> combinations = _configuration->getAllRealCombinations();
            for(size_t i = 0; i < combinations.size(); i++){
                ++keys;
                SimulationKey k;
                KnobsValues v = combinations.at(i);
                k.numCores = v[KNOB_VIRTUAL_CORES];
                k.frequency = v[KNOB_FREQUENCY];

                primaryPrediction = sel->getPrimaryPrediction(v);
                secondaryPrediction = sel->getSecondaryPrediction(v);

                switch(_p.contractType){
                    case CONTRACT_PERF_UTILIZATION:
                    case CONTRACT_PERF_BANDWIDTH:
                    case CONTRACT_PERF_COMPLETION_TIME:{
                        primaryValue = 1.0 / table[k].completionTime;
                        secondaryValue = table[k].wattsCores;
                        double bwError = abs((primaryValue - primaryPrediction) / primaryPrediction)*100.0;
                        double pwError = abs((secondaryValue - secondaryPrediction) / secondaryPrediction)*100.0;
                        mapeBandwidth += bwError;
                        mapePower += pwError;
                        res.performanceErrors.push_back(bwError);
                        res.powerErrors.push_back(pwError);
                    }break;
                    case CONTRACT_POWER_BUDGET:{
                        secondaryValue = 1.0 / table[k].completionTime;
                        primaryValue = table[k].wattsCores;
                        double bwError = abs((secondaryValue - secondaryPrediction) / secondaryPrediction)*100.0;
                        double pwError = abs((primaryValue - primaryPrediction) / primaryPrediction)*100.0;
                        mapeBandwidth += bwError;
                        mapePower += pwError;
                        res.performanceErrors.push_back(bwError);
                        res.powerErrors.push_back(pwError);
                    }break;
                    default:{
                        ;
                    }break;
                }
            }

            mapeBandwidth /= keys;
            mapePower /= keys;

            res.bandwidthAccuracy = 100.0 - mapeBandwidth;
            res.powerAccuracy = 100.0 - mapePower;
        }break;
        default:{
            ;
        }break;
    }
    return res;
#endif
    return SimulationResult();
}

}
