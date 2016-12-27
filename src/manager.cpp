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

#include "manager.hpp"

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
#include <sys/wait.h>

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

static Parameters& validate(Parameters& p){
    ParametersValidation apv = p.validate();
    if(apv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + std::to_string(apv));
    }
    return p;
}

Manager::Manager(Parameters adaptivityParameters):
        _terminated(false),
        _p(validate(adaptivityParameters)),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _counter(_p.mammut.getInstanceEnergy()->getCounter()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _samples(initSamples()),
        _variations(new MovingAverageExponential<double>(0.5)),
        _totalTasks(0),
        _remainingTasks(0),
        _deadline(0),
        _lastStoredSampleMs(0),
        _inhibited(false),
        _configuration(NULL),
        _selector(NULL),
        _pid(0)
{
    _cpufreq->removeTurboFrequencies();
}

Manager::~Manager(){
    ;
}

void Manager::run(){
    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        _remainingTasks = _p.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 + _p.requiredCompletionTime;
    }else if(_p.contractType == CONTRACT_PERF_MAX){
        // In this way all the configurations will be suboptimal and
        // we will select the less suboptimal one, i.e. the closest
        // to the requirement, i.e. the most performing one.
        _p.requiredBandwidth = std::numeric_limits<double>::max();
    }

    waitForStart();
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

    /** Wait for the disinhibition from global manager. **/
    while(_inhibited){;}

    if(_counter){
        _counter->reset();
    }
    //TODO RESET BANDWIDTHIN
    _lastStoredSampleMs = getMillisecondsTime();
    if(_p.observer){
        _p.observer->_startMonitoringMs = _lastStoredSampleMs;
    }

    /* Force the first calibration point. **/
    KnobsValues kv = decide();
    act(kv);

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
        if(!_inhibited){
            observe();
            updateRequiredBandwidth();
            logObservation();

            if(!persist()){
                DEBUG("Asking selector.");
                KnobsValues kv = decide();
                act(kv);
                _configuration->trigger();
                startSample = getMillisecondsTime();
            }
        }else{
            // If inhibited, we need to discard the previous samples
            // (they were taken with less/more applications running).
            _samples->reset();
        }
    }

    clean();
    ulong duration = getExecutionTime();
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
        ReconfigurationStats rs = _configuration->getReconfigurationStats();
        _p.observer->summaryStats(cs, rs, duration, _totalTasks);
    }
}

void Manager::terminate(){
    _terminated = true;
}

void Manager::inhibit(){
    _inhibited = true;
}

void Manager::disinhibit(){
    _inhibited = false;
}

void Manager::shrink(CalibrationShrink type){
    switch(type){
        case CALIBRATION_SHRINK_AGGREGATE:{
            if(_pid){
                _task->getProcessHandler(_pid)->move(_topology->getVirtualCores().back());
            }
        }break;
        case CALIBRATION_SHRINK_PAUSE:{
            shrinkPause();
        }break;
        default:{
            ;
        }break;
    }
}

void Manager::stretch(CalibrationShrink type){
    switch(type){
        case CALIBRATION_SHRINK_AGGREGATE:{
            if(_pid){
                std::vector<VirtualCore*> allowedCores = ((KnobMapping*) _configuration->getKnob(KNOB_TYPE_MAPPING))->getAllowedCores();
                _task->getProcessHandler(_pid)->move(allowedCores);
            }
        }break;
        case CALIBRATION_SHRINK_PAUSE:{
            stretchPause();
        }break;
        default:{
            ;
        }break;
    }
}

void Manager::updateModelsInterference(){
    if(_p.strategySelection != STRATEGY_SELECTION_LEARNING){
        throw std::runtime_error("updateModelsInterference can only be "
                "used on LEARNING selectors.");
    }
    // Temporarely disinhibit.
    _samples->reset();
    disinhibit();
    ((SelectorLearner*) _selector)->updateModelsInterference();
}

void Manager::waitModelsInterferenceUpdate(){
    if(_p.strategySelection != STRATEGY_SELECTION_LEARNING){
        throw std::runtime_error("waitModelsInterferenceUpdate can only be "
                        "used on LEARNING selectors.");
    }
    while(!((SelectorLearner*) _selector)->areModelsUpdated() && !_terminated){
        // If in the meanwhile the selector is waiting for calibration,
        // allow him to calibrate.
        if(_selector->isCalibrating()){
           _selector->allowCalibration();
        }
    }
    // Inhibit again.
    inhibit();
}

std::vector<VirtualCoreId> Manager::getUsedCores(){
    while(_selector->isCalibrating() || _selector->getTotalCalibrationTime() == 0){;}
    vector<VirtualCore*> vc = ((KnobMapping*) _configuration->getKnob(KNOB_TYPE_MAPPING))->getActiveVirtualCores();
    vector<VirtualCoreId> vcid;
    for(VirtualCore* v : vc){
        vcid.push_back(v->getVirtualCoreId());
    }
    return vcid;
}

void Manager::allowCores(std::vector<mammut::topology::VirtualCoreId> ids){
    vector<VirtualCore*> allowedVc;
    for(auto id : ids){
        allowedVc.push_back(_topology->getVirtualCore(id));
    }
    ((KnobVirtualCores*) _configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES))->changeMax(ids.size() - _configuration->getNumServiceNodes());
    ((KnobMapping*) _configuration->getKnob(KNOB_TYPE_MAPPING))->setAllowedCores(allowedVc);
}

void Manager::updateRequiredBandwidth() {
    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        double now = getMillisecondsTime();
        if(now / 1000.0 >= _deadline){
            _p.requiredBandwidth = numeric_limits<double>::max();
        }else{
            _p.requiredBandwidth = _remainingTasks / ((_deadline * 1000.0 - now) / 1000.0);
        }
    }
}

void Manager::setDomainToHighestFrequency(const Domain* domain){
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

double Manager::getPrimaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_MAX:{
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

double Manager::getSecondaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return sample.watts;
        }break;
        case CONTRACT_PERF_MAX:{
            // In this case we just ignore the power consumption.
            // By putting it to 0, we force the algorithm to select
            // the most performing configuration.
            return 0;
        }
        case CONTRACT_POWER_BUDGET:{
            return sample.bandwidth;
        }break;
        default:{
            return 0;
        }break;
    }
}

bool Manager::persist() const{
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

void Manager::lockKnobs() const{
    if(!_p.knobCoresEnabled){
        _configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->lockToMax();
    }
    if(!_p.knobMappingEnabled){
        _configuration->getKnob(KNOB_TYPE_MAPPING)->lock(MAPPING_TYPE_LINEAR);
    }
    if(!_p.knobFrequencyEnabled){
        _configuration->getKnob(KNOB_TYPE_FREQUENCY)->lockToMax();
    }
    if(!_p.knobHyperthreadingEnabled){
        _configuration->getKnob(KNOB_TYPE_HYPERTHREADING)->lockToMin();
    }
}

Selector* Manager::createSelector() const{
    if(_p.contractType == CONTRACT_NONE){
        return new SelectorFixed(_p, *_configuration, _samples);
    }else{
        switch(_p.strategySelection){
            case STRATEGY_SELECTION_ANALYTICAL:{
                return new SelectorAnalytical(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LEARNING:{
                return new SelectorLearner(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LIMARTINEZ:{
                return new SelectorLiMartinez(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_MISHRA:{
                return new SelectorMishra(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_FULLSEARCH:{
                return new SelectorFullSearch(_p, *_configuration, _samples);
            }break;
            default:{
                throw std::runtime_error("Selector not yet implemented.");
            }break;
        }
    }
    return NULL;
}

Smoother<MonitoredSample>* Manager::initSamples() const{
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

void Manager::updateTasksCount(orlog::ApplicationSample& sample){
    double newTasks = sample.tasksCount;
    if(_p.synchronousWorkers){
        // When we have synchronous workers we need to count the iterations,
        // not the real tasks (indeed in this case each worker will receive
        // the same amount of tasks, e.g. in canneal) since they are sent in
        // broadcast.
        newTasks /= _configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealValue();
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

void Manager::observe(){
    MonitoredSample sample;
    orlog::ApplicationSample ws;
    Joules joules = 0.0;

    askSample();
    getSample(ws);
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
        sample.bandwidth /= _configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealValue();
    }

    _samples->add(sample);
    _variations->add(getPrimaryValue(_samples->coefficientVariation()));

    DEBUGB(samplesFile << *_samples << "\n");
}

KnobsValues Manager::decide(){
    if(!_configuration->knobsChangeNeeded()){
        return _configuration->getRealValues();
    }

    if(_totalTasks){
        // We need to update the input bandwidth only if we already processed
        // some tasks.  By doing this check, we avoid updating bandwidth for
        // the first forced reconfiguration.
        _selector->updateBandwidthIn();
    }
    return _selector->getNextKnobsValues(_totalTasks);
}

void Manager::act(KnobsValues kv, bool force){
    if(force || !_configuration->equal(kv)){
        _configuration->setValues(kv);
        manageConfigurationChange();

        /****************** Clean state ******************/
        _samples->reset();
        _variations->reset();
        DEBUG("Resetting sample.");
        // _lastStoredSampleMs = getMillisecondsTime();
        // Discarding a sample since it may be between 2 different configurations.
        orlog::ApplicationSample ws;
        askSample();
        getSample(ws);
        updateTasksCount(ws);

        Joules joules = getAndResetJoules();
        if(_p.observer){
            _p.observer->addJoules(joules);
        }
        _lastStoredSampleMs = getMillisecondsTime();
    }
}

Joules Manager::getAndResetJoules(){
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

void Manager::logObservation(){
    if(_p.observer){
        const KnobMapping* kMapping = dynamic_cast<const KnobMapping*>(_configuration->getKnob(KNOB_TYPE_MAPPING));
        MonitoredSample ms = _samples->average();
        _p.observer->observe(_lastStoredSampleMs,
                             _configuration->getRealValue(KNOB_TYPE_VIRTUAL_CORES),
                             _configuration->getRealValue(KNOB_TYPE_FREQUENCY),
                             kMapping->getActiveVirtualCores(),
                             _samples->getLastSample().bandwidth,
                             ms.bandwidth,
                             _samples->coefficientVariation().bandwidth,
                             ms.latency,
                             ms.utilisation,
                             _samples->getLastSample().watts,
                             ms.watts);
    }
}

ManagerExternal::ManagerExternal(const std::string& orlogChannel,
                                 Parameters adaptivityParameters):
        Manager(adaptivityParameters), _monitor(orlogChannel){
    Manager::_configuration = new ConfigurationExternal(_p);
    lockKnobs();
    _configuration->createAllRealCombinations();
    _selector = createSelector();
    // For external application we do not care if synchronous of not (we count iterations).   
    _p.synchronousWorkers = false;
}

ManagerExternal::ManagerExternal(nn::socket& orlogSocket,
                                 int chid,
                                 Parameters adaptivityParameters):
            Manager(adaptivityParameters), _monitor(orlogSocket, chid){
    Manager::_configuration = new ConfigurationExternal(_p);
    lockKnobs();
    _configuration->createAllRealCombinations();
    _selector = createSelector();
    // For external application we do not care if synchronous of not (we count iterations).
    _p.synchronousWorkers = false;
}

ManagerExternal::~ManagerExternal(){
    if(Manager::_configuration){
        delete Manager::_configuration;
    }
    if(_selector){
        delete _selector;
    }
}

void ManagerExternal::waitForStart(){
    Manager::_pid = _monitor.waitStart();
    ((KnobMappingExternal*)_configuration->getKnob(KNOB_TYPE_MAPPING))->setPid(_pid);
}

void ManagerExternal::askSample(){
    return;
}
void ManagerExternal::getSample(orlog::ApplicationSample& sample){
    if(!_monitor.getSample(sample)){
        _terminated = true;
    }
}

void ManagerExternal::manageConfigurationChange(){;}

void ManagerExternal::clean(){;}

ulong ManagerExternal::getExecutionTime(){
    return _monitor.getExecutionTime();
}

void ManagerExternal::shrinkPause(){
    kill(_pid, SIGSTOP);
}

void ManagerExternal::stretchPause(){
    kill(_pid, SIGCONT);
}

ManagerBlackBox::ManagerBlackBox(pid_t pid, Parameters adaptivityParameters):
            Manager(adaptivityParameters), _process(adaptivityParameters.mammut.getInstanceTask()->getProcessHandler(pid)),
            _startTime(getMillisecondsTime()), _lastTime(_startTime){
        Manager::_pid = pid;
        Manager::_configuration = new ConfigurationExternal(_p);
        lockKnobs();
        _configuration->createAllRealCombinations();
        _selector = createSelector();
        // For external application we do not care if synchronous of not (we count iterations).
        _p.synchronousWorkers = false;
        // Check supported contracts.
        if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME ||
           _p.contractType == CONTRACT_PERF_BANDWIDTH ||
           _p.contractType == CONTRACT_PERF_UTILIZATION){
            throw std::runtime_error("ManagerBlackBox. Unsupported contract.");
        }
}

ManagerBlackBox::~ManagerBlackBox(){
    if(Manager::_configuration){
        delete Manager::_configuration;
    }
    if(_selector){
        delete _selector;
    }
}

pid_t ManagerBlackBox::getPid() const{
    return _pid;
}

void ManagerBlackBox::waitForStart(){
    // We know for sure that when the Manager is created the process
    // already started.
    ((KnobMappingExternal*)_configuration->getKnob(KNOB_TYPE_MAPPING))->setProcessHandler(_process);
}
void ManagerBlackBox::askSample(){
    // We do not need to ask for black box.
    ;
}

static bool isRunning(pid_t pid) {
    return !kill(pid, 0);
}

void ManagerBlackBox::getSample(orlog::ApplicationSample& sample){
    double now = getMillisecondsTime();
    double cycles;
    if(!isRunning(_process->getId())){
        _terminated = true;
        return;
    }
    _process->getAndResetCycles(cycles);
    sample.bandwidthTotal = cycles / (now - _lastTime);
    sample.latency = -1; // Not used.
    sample.loadPercentage = 100.0; // We do not know what's the input bandwidth.
    sample.tasksCount = 0; // We do not know how many iterations have been performed.
    _lastTime = now;
}

void ManagerBlackBox::manageConfigurationChange(){;}

void ManagerBlackBox::clean(){;}

ulong ManagerBlackBox::getExecutionTime(){
    return getMillisecondsTime() - _startTime;
}

void ManagerBlackBox::shrinkPause(){
    pid_t pid = _process->getId();
    kill(pid, SIGSTOP);
}

void ManagerBlackBox::stretchPause(){
    pid_t pid = _process->getId();
    kill(pid, SIGCONT);
}

}
