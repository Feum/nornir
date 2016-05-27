/*
 * manager.cpp
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
#include "./manager.hpp"
#include "./parameters.hpp"
#include "./predictors.hpp"
#include "./node.hpp"
#include "./utils.hpp"

#include "external/Mammut/mammut/module.hpp"
#include "external/Mammut/mammut/utils.hpp"
#include "external/Mammut/mammut/mammut.hpp"

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
void ManagerFarm<lb_t, gt_t>::setDomainToHighestFrequency(const Domain* domain){
    if(!domain->setGovernor(GOVERNOR_PERFORMANCE)){
        if(!domain->setGovernor(GOVERNOR_USERSPACE) ||
           !domain->setHighestFrequencyUserspace()){
            throw runtime_error("AdaptivityManagerFarm: Fatal error while "
                                "setting highest frequency for sensitive "
                                "emitter/collector. Try to run it without "
                                "sensitivity parameters.");
        }
    }
}

template <typename lb_t, typename gt_t>
double ManagerFarm<lb_t, gt_t>::getPrimaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return sample.bandwidth;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return sample.watts;
        }break;
        default:{
            return 0;
        }break;
    }
}

template <typename lb_t, typename gt_t>
double ManagerFarm<lb_t, gt_t>::getSecondaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return sample.watts;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return sample.bandwidth;
        }break;
        default:{
            return 0;
        }break;
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::changeKnobs(){
    if(!_configuration.knobsChangeNeeded()){
        return;
    }

    KnobsValues values = _selector->getNextKnobsValues(_totalTasks);
    if(!_configuration.equal(values)){
        _configuration.setValues(values);

        const KnobWorkers* knobWorkers = dynamic_cast<const KnobWorkers*>(_configuration.getKnob(KNOB_TYPE_WORKERS));
        std::vector<AdaptiveNode*> newWorkers = knobWorkers->getActiveWorkers();
        WorkerSample ws;

        if(_activeWorkers.size() != newWorkers.size()){
            /** 
             * Since I stopped the workers after I asked for a sample, there
             * may still be tasks that have been processed but I did not count.
             * For this reason, I get them.
             * I do not need to ask since the node put it in the Q when it 
             * terminated.
             */
            DEBUG("Getting spurious..");
            getWorkersSamples(ws);
            updateTasksCount(ws);
            DEBUG("Spurious got.");
        }

        _activeWorkers = newWorkers;

        /****************** Clean state ******************/
        _samples->reset();
        _variations->reset();
        Joules joules = getAndResetJoules();
        if(_p.observer){
            _p.observer->addJoules(joules);
        }
        DEBUG("Resetting sample.");
        _lastStoredSampleMs = getMillisecondsTime();
        askForWorkersSamples();
        getWorkersSamples(ws);
        updateTasksCount(ws);
        //resetSample();
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::observe(){
    if(_p.observer){
        const KnobMapping* kMapping = dynamic_cast<const KnobMapping*>(_configuration.getKnob(KNOB_TYPE_MAPPING));
        MonitoredSample ms = _samples->average();
        _p.observer->observe(_lastStoredSampleMs,
                             _configuration.getRealValue(KNOB_TYPE_WORKERS),
                             _configuration.getRealValue(KNOB_TYPE_FREQUENCY),
                             kMapping->getEmitterVirtualCore(),
                             kMapping->getWorkersVirtualCore(),
                             kMapping->getCollectorVirtualCore(),
                             _samples->getLastSample().bandwidth,
                             ms.bandwidth,
                             _samples->coefficientVariation().bandwidth,
                             ms.latency,
                             ms.utilisation,
                             _samples->getLastSample().watts,
                             ms.watts);
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::askForWorkersSamples(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->askForSample();
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::getWorkersSamples(WorkerSample& sample){
    AdaptiveNode* w;
    uint numActiveWorkers = _activeWorkers.size();
    sample = WorkerSample();

    for(size_t i = 0; i < numActiveWorkers; i++){
        WorkerSample tmp;
        w = _activeWorkers.at(i);
        w->getSampleResponse(tmp, _samples->average().latency);
        sample += tmp;
    }
    sample.loadPercentage /= numActiveWorkers;
    sample.latency /= numActiveWorkers;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::resetSample(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->resetSample();
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::updateTasksCount(WorkerSample& sample){
    double newTasks = sample.tasksCount;
    if(_p.synchronousWorkers){
        // When we have synchronous workers we need to count the iterations,
        // not the real tasks (indeed in this case each worker will receive
        // the same amount of tasks, e.g. in canneal) since they are sent in
        // broadcast.
        newTasks /= _activeWorkers.size();
    }
    _totalTasks += newTasks;
    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        if(_remainingTasks > newTasks){
            _remainingTasks -= newTasks;
        }else{
            _remainingTasks = 0;
        }
    }
}

template <typename lb_t, typename gt_t>
Joules ManagerFarm<lb_t, gt_t>::getAndResetJoules(){
    Joules joules = 0.0;
    if(_counter){
        switch(_counter->getType()){
            case COUNTER_CPUS:{
                joules = ((CounterCpus*) _counter)->getJoulesCoresAll();
            }break;
            default:{
                joules = _counter->getJoules();
            }break;
        }
        _counter->reset();
    }
    return joules;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::storeNewSample(){
    MonitoredSample sample;
    WorkerSample ws;
    Joules joules = 0.0;

    askForWorkersSamples();
    getWorkersSamples(ws);
    updateTasksCount(ws);

    joules = getAndResetJoules();
    if(_p.observer){
        _p.observer->addJoules(joules);
    }

    double now = getMillisecondsTime();
    double durationSecs = (now - _lastStoredSampleMs) / 1000.0;
    _lastStoredSampleMs = now;

    sample.watts = joules / durationSecs;
    // ATTENTION: Bandwidth is not the number of task since the
    //            last observation but the number of expected
    //            tasks that will be processed in 1 second.
    //            For this reason, if we sum all the bandwidths in
    //            the result observation file, we may have an higher
    //            number than the number of tasks.
    sample.bandwidth = ws.bandwidthTotal;
    sample.utilisation = ws.loadPercentage;
    sample.latency = ws.latency;

    if(_p.synchronousWorkers){
        // When we have synchronous workers we need to divide
        // for the number of workers since we do it for the totalTasks
        // count. When this flag is set we count iterations, not real
        // tasks.
        sample.bandwidth /= _activeWorkers.size();
    }

    _samples->add(sample);
    _variations->add(getPrimaryValue(_samples->coefficientVariation()));

    DEBUGB(samplesFile << *_samples << "\n");
}

template <typename lb_t, typename gt_t>
bool ManagerFarm<lb_t, gt_t>::persist() const{
    bool r = false;
    switch(_p.strategyPersistence){
        case STRATEGY_PERSISTENCE_SAMPLES:{
            r = _samples->size() < _p.persistenceValue;
        }break;
        case STRATEGY_PERSISTENCE_VARIATION:{
            const MonitoredSample& variation = _samples->coefficientVariation();
            double primaryVariation =  getPrimaryValue(variation);
            double secondaryVariation =  getSecondaryValue(variation);
            r = _samples->size() < 1 ||
                primaryVariation > _p.persistenceValue ||
                secondaryVariation > _p.persistenceValue;
#if 0
            double primaryVariation = _variations->coefficientVariation();
            std::cout << "Variation size: " << _variations->size() << " PrimaryVariation: " << primaryVariation << std::endl;
            r = _variations->size() < 2 || primaryVariation > _p.persistenceValue;
#endif
        }break;
    }
    return r;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initSelector(){
    if(_p.contractType == CONTRACT_NONE){
        _selector = new SelectorFixed(_p, _configuration, _samples);
    }else{
        switch(_p.strategySelection){
            case STRATEGY_SELECTION_ANALYTICAL:{
                _selector = new SelectorAnalytical(_p, _configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LEARNING:{
                _selector = new SelectorLearner(_p, _configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LIMARTINEZ:{
                _selector = new SelectorLiMartinez(_p, _configuration, _samples);
            }break;
            case STRATEGY_SELECTION_MISHRA:{
                _selector = new SelectorMishra(_p, _configuration, _samples);
            }break;
            case STRATEGY_SELECTION_FULLSEARCH:{
                _selector = new SelectorFullSearch(_p, _configuration, _samples);
            }break;
            default:{
                throw std::runtime_error("Selector not yet implemented.");
            }break;
        }
    }
}

static Parameters& validate(Parameters& p){
    ParametersValidation apv = p.validate();
    if(apv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + std::to_string(apv));
    }
    return p;
}

static std::vector<AdaptiveNode*> convertWorkers(svector<ff_node*> w){
    std::vector<AdaptiveNode*> r;
    for(size_t i = 0; i < w.size(); i++){
        r.push_back(dynamic_cast<AdaptiveNode*>(w[i]));
    }
    return r;
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::ManagerFarm(ff_farm<lb_t, gt_t>* farm, Parameters parameters):
        _farm(farm),
        _terminated(false),
        _p(validate(parameters)),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _counter(_p.mammut.getInstanceEnergy()->getCounter()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _emitter(dynamic_cast<AdaptiveNode*>(_farm->getEmitter())),
        _collector(dynamic_cast<AdaptiveNode*>(_farm->getCollector())),
        _activeWorkers(convertWorkers(_farm->getWorkers())),
        _samples(initSamples()),
        _variations(new MovingAverageExponential<double>(0.5)),
        _configuration(_p, _emitter, _collector, _farm->getgt(), _activeWorkers,
                       _samples, &_terminated),
        _totalTasks(0),
        _remainingTasks(0),
        _deadline(0),
        _lastStoredSampleMs(0),
        _selector(NULL){
    DEBUGB(samplesFile.open("samples.csv"));
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::~ManagerFarm(){
    delete _samples;
    delete _variations;
    if(_selector){
        delete _selector;
    }
    DEBUGB(samplesFile.close());
}

template <typename lb_t, typename gt_t>
Smoother<MonitoredSample>* ManagerFarm<lb_t, gt_t>::initSamples() const{
    switch(_p.strategySmoothing){
        case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
            return new MovingAverageSimple<MonitoredSample>(_p.smoothingFactor);
        }break;
        case STRATEGY_SMOOTHING_EXPONENTIAL:{
            return new MovingAverageExponential<MonitoredSample>(_p.smoothingFactor);
        }break;
        default:{
            return NULL;
        }
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
void ManagerFarm<lb_t, gt_t>::cleanNodes() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->clean();
    }
    if (_emitter) {
        _emitter->clean();
    }
    if (_collector) {
        _collector->clean();
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::updateRequiredBandwidth() {
    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        double now = getMillisecondsTime();
        if(now / 1000.0 >= _deadline){
            _p.requiredBandwidth = numeric_limits<double>::max();
        }else{
            _p.requiredBandwidth = _remainingTasks / ((_deadline * 1000.0 - now) / 1000.0);
        }
    }
}

#define PAR_BEGIN_ENV "__PAR_BEGIN"

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
SimulationResult ManagerFarm<lb_t, gt_t>::simulate(string configurationData, volatile bool* terminate, size_t maxIterations){
    vector<string> lines = readFile(configurationData);
    map<SimulationKey, SimulationData, SimulationKeyCompare> table;
    KnobsValues lastConfigurationValues = _configuration.getRealValues();
    // Starts from 1 to skip the header line.
    for(size_t i = 1; i < lines.size(); i++){
        SimulationKey key;
        SimulationData data;
        vector<string> fields = split(lines.at(i), '\t');
        key.numCores = atof(fields[0].c_str());
        key.frequency = atof(fields[1].c_str());

        data.completionTime = atof(fields[2].c_str());
        data.wattsCpu = atof(fields[3].c_str());
        data.wattsCores = atof(fields[4].c_str());
        data.wattsDram = atof(fields[5].c_str());

        table.insert(std::pair<SimulationKey, SimulationData>(key, data));
    }

    initSelector();

    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        _remainingTasks = _p.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 + _p.requiredCompletionTime;
    }

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
#if 0
    /** Creates the parallel section begin file. **/
    char* default_in_roi = (char*) malloc(sizeof(char)*256);
    default_in_roi[0] = '\0';
    default_in_roi = strcat(default_in_roi, getenv("HOME"));
    default_in_roi = strcat(default_in_roi, "/roi_in");
    setenv(PAR_BEGIN_ENV, default_in_roi, 0);
    free(default_in_roi);
    FILE* in_roi = fopen(getenv(PAR_BEGIN_ENV), "w");
    fclose(in_roi);
#endif

    if(_counter){
        _counter->reset();
    }
    _lastStoredSampleMs = getMillisecondsTime();
    if(_p.observer){
        _p.observer->_startMonitoringMs = _lastStoredSampleMs;
    }

    /* Force the first calibration point. **/
    changeKnobs();

    _samples->reset();
    _variations->reset();

    SimulationResult res;
    KnobsValues values = _configuration.getRealValues();
    SimulationKey key;
    SimulationData data;
    key.numCores = values[KNOB_TYPE_WORKERS];
    key.frequency = values[KNOB_TYPE_FREQUENCY];
    if(table.find(key) != table.end()){
        data = table[key];
    }else{
        throw std::runtime_error("Impossible to find value for key.");
    }

    MonitoredSample sample;
    sample.watts = data.wattsCores;
    sample.bandwidth = 1.0 / data.completionTime;
    sample.utilisation = 100.0; //TODO Not yet implemented.
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

    while(!_configuration.equal(lastConfigurationValues) && (!maxIterations || steps <= maxIterations)){
        ++steps;
        lastConfigurationValues = _configuration.getRealValues();
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
        storeNewSample();
        _samples->reset();
        _variations->reset();


        MonitoredSample sample;
        sample.watts = data.wattsCores;
        sample.bandwidth = 1.0 / data.completionTime;
        sample.utilisation = 100.0; //TODO Not yet implemented.
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

        observe();

        if(!persist()){
            DEBUG("Asking selector.");
            changeKnobs();

            values = _configuration.getRealValues();
            key.numCores = values[KNOB_TYPE_WORKERS];
            key.frequency = values[KNOB_TYPE_FREQUENCY];
            if(table.find(key) != table.end()){
                data = table[key];
            }else{
                throw std::runtime_error("Impossible to find value for key.");
            }

            _configuration.trigger();
            startSample = getMillisecondsTime();
        }
    }

#if 0
    unlink(getenv(PAR_BEGIN_ENV));
#endif

    *terminate = true;
    DEBUG("Terminating...wait freezing.");
    _farm->wait_freezing();
    _farm->wait();
    DEBUG("Terminated.");

    cleanNodes();

    // We do -1 because the last step was the optimal configuration.
    res.numSteps = steps - 1;
    res.foundConfiguration = _configuration.getRealValues();

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
            std::vector<KnobsValues> combinations = _configuration.getAllRealCombinations();
            for(size_t i = 0; i < combinations.size(); i++){
                ++keys;
                SimulationKey k;
                KnobsValues v = combinations.at(i);
                k.numCores = v[KNOB_TYPE_WORKERS];
                k.frequency = v[KNOB_TYPE_FREQUENCY];

                primaryPrediction = sel->getPrimaryPrediction(v);
                secondaryPrediction = sel->getSecondaryPrediction(v);

                switch(_p.contractType){
                    case CONTRACT_PERF_UTILIZATION:
                    case CONTRACT_PERF_BANDWIDTH:
                    case CONTRACT_PERF_COMPLETION_TIME:{
                        primaryValue = 1.0 / table[k].completionTime;
                        secondaryValue = table[k].wattsCores;
                        mapeBandwidth += abs((primaryValue - primaryPrediction) / primaryPrediction)*100.0;
                        mapePower += abs((secondaryValue - secondaryPrediction) / secondaryPrediction)*100.0;
                    }break;
                    case CONTRACT_POWER_BUDGET:{
                        secondaryValue = 1.0 / table[k].completionTime;
                        primaryValue = table[k].wattsCores;
                        mapePower += abs((secondaryValue - secondaryPrediction) / secondaryPrediction)*100.0;
                        mapePower += abs((primaryValue - primaryPrediction) / primaryPrediction)*100.0;
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
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::run(){
    initSelector();

    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        _remainingTasks = _p.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 + _p.requiredCompletionTime;
    }

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
#if 0
    /** Creates the parallel section begin file. **/
    char* default_in_roi = (char*) malloc(sizeof(char)*256);
    default_in_roi[0] = '\0';
    default_in_roi = strcat(default_in_roi, getenv("HOME"));
    default_in_roi = strcat(default_in_roi, "/roi_in");
    setenv(PAR_BEGIN_ENV, default_in_roi, 0);
    free(default_in_roi);
    FILE* in_roi = fopen(getenv(PAR_BEGIN_ENV), "w");
    fclose(in_roi);
#endif

    if(_counter){
        _counter->reset();
    }
    _lastStoredSampleMs = getMillisecondsTime();
    if(_p.observer){
        _p.observer->_startMonitoringMs = _lastStoredSampleMs;
    }

    /* Force the first calibration point. **/
    changeKnobs();

    ThreadHandler* thisThread = _task->getProcessHandler(getpid())->getThreadHandler(gettid());
    thisThread->move(MANAGER_VIRTUAL_CORE);

    double microsecsSleep = 0;
    double startSample = getMillisecondsTime();
    double overheadMs = 0;

    uint samplingInterval;
    uint steadySamples = 0;

    while(!_terminated){
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
        storeNewSample();
        DEBUG("New sample stored.");

        updateRequiredBandwidth();

        observe();

        if(!persist()){
            DEBUG("Asking selector.");
            changeKnobs();
            _configuration.trigger();
            startSample = getMillisecondsTime();
        }
    }

    DEBUG("Terminating...wait freezing.");
    _farm->wait_freezing();
    _farm->wait();
    DEBUG("Terminated.");

    double duration = _farm->ffTime();
#if 0
    unlink(getenv(PAR_BEGIN_ENV));
#endif
    if(_p.observer){
        vector<CalibrationStats> cs;
        if(_selector){
            _selector->stopCalibration(_totalTasks);
            cs = _selector->getCalibrationsStats();
            _p.observer->calibrationStats(cs, duration, _totalTasks);
        }
        ReconfigurationStats rs = _configuration.getReconfigurationStats();
        _p.observer->summaryStats(cs, rs, duration, _totalTasks);
    }

    cleanNodes();
}

}
