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
        throw runtime_error("Invalid parameters: " + std::to_string(apv));
    }
    return p;
}

Manager::Manager(Parameters nornirParameters):
        _terminated(false),
        _p(validate(nornirParameters)),
        _counter(NULL),
        _task(NULL),
        _topology(NULL),
        _samples(initSamples()),
        _variations(new MovingAverageExponential<double>(0.5)),
        _totalTasks(0),
        _remainingTasks(0),
        _deadline(0),
        _lastStoredSampleMs(0),
        _inhibited(false),
        _configuration(NULL),
        _selector(NULL),
        _pid(0),
        _toSimulate(false)
{
    DEBUG("Initializing manager.");
    for(LoggerType lt : _p.loggersTypes){
        switch(lt){
        case LOGGER_FILE:{
            _p.loggers.push_back(new LoggerFile());
        }break;
        case LOGGER_GRAPHITE:{
            // TODO passare comeparametro url e porta
            _p.loggers.push_back(new LoggerGraphite("178.62.14.35", 2003));
        }break;
        default:{
            throw std::runtime_error("Unknown logger type.");
        }
        }
    }
    DEBUGB(samplesFile.open("samples.csv"));
}

Manager::~Manager(){
    for(auto logger : _p.loggers){
        delete logger; //TODO: NO, I may also delete loggers created by the user.
    }
    _p.loggers.clear();
    DEBUGB(samplesFile.close());
}

void Manager::run(){
    _p.mammut.getInstanceCpuFreq()->removeTurboFrequencies();
    DEBUG("Turbo frequencies removed.");
    if(!_toSimulate){
        _counter = _p.mammut.getInstanceEnergy()->getCounter();
        _task = _p.mammut.getInstanceTask();
        _topology = _p.mammut.getInstanceTopology();
    }
    DEBUG("Mammut handlers created.");

    if(isPrimaryRequirement(_p.requirements.executionTime)){
        _remainingTasks = _p.requirements.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 + _p.requirements.executionTime;
    }

    waitForStart();
    lockKnobs();
    _configuration->createAllRealCombinations();
    _selector = createSelector();
    for(auto logger : _p.loggers){
        logger->setStartTimestamp();
    }
    DEBUG("Application started.");

    /** Wait for the disinhibition from global manager. **/
    while(_inhibited){;}

    // Reset joules counter
    getAndResetJoules();
    //TODO RESET BANDWIDTHIN
    _lastStoredSampleMs = getMillisecondsTime();

    /* Force the first calibration point. **/
    KnobsValues kv = decide();
    act(kv);

    if(!_toSimulate){
        ThreadHandler* thisThread = _task->getProcessHandler(getpid())->getThreadHandler(gettid());
        thisThread->move(NORNIR_MANAGER_VIRTUAL_CORE);
    }

    double startSample = getMillisecondsTime();
    uint samplingInterval, steadySamples = 0;

    while(!_terminated){
        double overheadMs = getMillisecondsTime() - startSample;
        if(_selector && _selector->isCalibrating()){
            samplingInterval = _p.samplingIntervalCalibration;
            steadySamples = 0;
        }else if(steadySamples < _p.steadyThreshold){
            samplingInterval = _p.samplingIntervalCalibration;
            steadySamples++;
        }else{
            samplingInterval = _p.samplingIntervalSteady;
        }
        double microsecsSleep = ((double)samplingInterval - overheadMs)*
                                 (double)MAMMUT_MICROSECS_IN_MILLISEC;
        if(microsecsSleep < 0){
            microsecsSleep = 0;
        }else{
            usleep(microsecsSleep);
        }

        startSample = getMillisecondsTime();
        if(!_inhibited){
            observe();
            updateRequiredThroughput();

            if(!persist()){
                DEBUG("Asking selector.");
                KnobsValues kv = decide();
                act(kv);
                startSample = getMillisecondsTime();
            }
        }else{
            // If inhibited, we need to discard the previous samples
            // (they were taken with less/more applications running).
            _samples->reset();
        }
    }

    terminationManagement(); // e.g. to collect final tasks count from riff
    ulong duration = getExecutionTime();

    for(auto logger : _p.loggers){
        logger->logSummary(*_configuration, _selector, duration, _totalTasks);
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
                std::vector<VirtualCore*> allowedCores = dynamic_cast<KnobMapping*>(_configuration->getKnob(KNOB_MAPPING))->getAllowedCores();
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
    dynamic_cast<SelectorLearner*>(_selector)->updateModelsInterference();
}

void Manager::waitModelsInterferenceUpdate(){
    if(_p.strategySelection != STRATEGY_SELECTION_LEARNING){
        throw std::runtime_error("waitModelsInterferenceUpdate can only be "
                        "used on LEARNING selectors.");
    }
    while(!dynamic_cast<SelectorLearner*>(_selector)->areModelsUpdated() && !_terminated){
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
    vector<VirtualCore*> vc = dynamic_cast<KnobMapping*>(_configuration->getKnob(KNOB_MAPPING))->getActiveVirtualCores();
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
    dynamic_cast<KnobVirtualCores*>(_configuration->getKnob(KNOB_VIRTUAL_CORES))->changeMax(ids.size() - _configuration->getNumServiceNodes());
    dynamic_cast<KnobMapping*>(_configuration->getKnob(KNOB_MAPPING))->setAllowedCores(allowedVc);
}

void Manager::setSimulationParameters(std::string samplesFileName){
    _toSimulate = true;
    std::ifstream file(samplesFileName);
    MonitoredSample sample;

    while(file >> sample){
        _simulationSamples.push_back(sample);
    }
}

void Manager::postConfigurationManagement(){;}

void Manager::terminationManagement(){;}

void Manager::updateRequiredThroughput() {
    if(isPrimaryRequirement(_p.requirements.executionTime)){
        double now = getMillisecondsTime();
        if(now / 1000.0 >= _deadline){
            _p.requirements.throughput = numeric_limits<double>::max();
        }else{
            _p.requirements.throughput = _remainingTasks / ((_deadline * 1000.0 - now) / 1000.0);
        }
    }
}

void Manager::setDomainToHighestFrequency(const Domain* domain){
    if(!domain->setGovernor(GOVERNOR_PERFORMANCE)){
        if(!domain->setGovernor(GOVERNOR_USERSPACE) ||
           !domain->setHighestFrequencyUserspace()){
            throw runtime_error("Manager: Fatal error when trying to set"
                                "domain to maximum frequency.");
        }
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
            r = _samples->size() < 1 ||
                variation.throughput > _p.persistenceValue ||
                variation.latency > _p.persistenceValue ||
                variation.watts > _p.persistenceValue;
        }break;
    }
    return r;
}

void Manager::lockKnobs() const{
    if(!_p.knobCoresEnabled){
        _configuration->getKnob(KNOB_VIRTUAL_CORES)->lockToMax();
    }
    if(!_p.knobMappingEnabled){
        _configuration->getKnob(KNOB_MAPPING)->lock(MAPPING_TYPE_LINEAR);
    }
    if(!_p.knobFrequencyEnabled){
        _configuration->getKnob(KNOB_FREQUENCY)->lockToMax();
    }
    if(!_p.knobHyperthreadingEnabled){
        _configuration->getKnob(KNOB_HYPERTHREADING)->lock(_p.knobHyperthreadingFixedValue);
    }
    if(!_p.knobClkModEmulatedEnabled){
        _configuration->getKnob(KNOB_CLKMOD_EMULATED)->lockToMax();
    }
}

Selector* Manager::createSelector() const{
    if(!_p.requirements.anySpecified() &&
        _p.strategySelection != STRATEGY_SELECTION_MANUAL_CLI &&
        _p.strategySelection != STRATEGY_SELECTION_MANUAL_WEB){
        // We use fixed selector if there were no requirements specified and
        // if the selection strategy is different from manual (indeed, for
        // manual selection there is no need to specify requirements since
        // the configuration is manually selected by some external entity).
        return new SelectorFixed(_p, *_configuration, _samples);
    }else{
        switch(_p.strategySelection){
            case STRATEGY_SELECTION_MANUAL_CLI:{
                return new SelectorManualCli(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_MANUAL_WEB:{
                return new SelectorManualWeb(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_ANALYTICAL:{
                return new SelectorAnalytical(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LEARNING:{
                return new SelectorLearner(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LIMARTINEZ:{
                return new SelectorLiMartinez(_p, *_configuration, _samples);
            }break;
            case STRATEGY_SELECTION_LEO:{
                return new SelectorLeo(_p, *_configuration, _samples);
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

void Manager::updateTasksCount(MonitoredSample &sample){
    _totalTasks += sample.numTasks;
    if(isPrimaryRequirement(_p.requirements.executionTime)){
        if(_remainingTasks > sample.numTasks){
            _remainingTasks -= sample.numTasks;
        }else{
            _remainingTasks = 0;
        }
    }
}

void Manager::observe(){
    Joules joules = 0.0;
    MonitoredSample sample;
    bool store = true;
    if(_toSimulate){
        if(_simulationSamples.empty()){
            return;
        }else{
            sample = _simulationSamples.front();
            _simulationSamples.pop_front();
            if(_simulationSamples.empty()){
                _terminated = true;
            }
        }
    }else{
        sample = getSample();
        if(_terminated){
            // When checking for sample, we may find out
            // that the application is terminated and we may
            // not have received an actual sample. In that
            // case, we do not store the sample (which is the last one).
            store = false;
        }else{
            double now = getMillisecondsTime();
            joules = getAndResetJoules();
            double durationSecs = (now - _lastStoredSampleMs) / 1000.0;
            _lastStoredSampleMs = now;

            // Add watts to the sample
            sample.watts = joules / durationSecs;

            if(_p.synchronousWorkers){
                // When we have synchronous workers we need to divide
                // for the number of workers since we do it for the totalTasks
                // count. When this flag is set we count iterations, not real
                // tasks.
                sample.throughput /= _configuration->getKnob(KNOB_VIRTUAL_CORES)->getRealValue();
                // When we have synchronous workers we need to count the iterations,
                // not the real tasks (indeed in this case each worker will receive
                // the same amount of tasks, e.g. in canneal) since they are sent in
                // broadcast.
                sample.numTasks /= _configuration->getKnob(KNOB_VIRTUAL_CORES)->getRealValue();
            }
        }
    }

    if(store){
        updateTasksCount(sample);
        _samples->add(sample);
        _variations->add(_samples->coefficientVariation().throughput);

        DEBUGB(samplesFile << sample << "\n");
        logObservation();
    }
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
        _selector->updateTotalTasks(_totalTasks);
    }
    return _selector->getNextKnobsValues();
}

void Manager::act(KnobsValues kv, bool force){
    if(force || !_configuration->equal(kv)){
        _configuration->setValues(kv);
        postConfigurationManagement();

        /****************** Clean state ******************/
        _samples->reset();
        _variations->reset();
        DEBUG("Resetting sample.");

        // We need to explicitely check that is not terminated.
        // Indeed, termination may have been detected by observe(),
        // but since we exit the control loop only after the act(),
        // at this point we may be still in the control loop
        // with _terminated = true.
        if(!_toSimulate && !_terminated){
            if(_p.cooldownPeriod){
                usleep(_p.cooldownPeriod * 1000);
            }
            // Don't store this sample since it may be inbetween 2
            // different configurations.
            MonitoredSample sample = clearStoredSample();
            updateTasksCount(sample);
            _lastStoredSampleMs = getMillisecondsTime();
            getAndResetJoules(); // Reset joules
        }
        _configuration->trigger();
    }
}

Joules Manager::getAndResetJoules(){
    Joules joules = 0.0;
    if(!_toSimulate){
        if(_counter){
            joules = _counter->getJoules();
            _counter->reset();
        }
    }
    return joules;
}

void Manager::logObservation(){
    for(auto logger : _p.loggers){
        logger->log(_selector->isCalibrating(), *_configuration, *_samples, _p.requirements);
    }
}

ManagerInstrumented::ManagerInstrumented(const std::string& riffChannel,
                                 Parameters nornirParameters):
        Manager(nornirParameters), _monitor(riffChannel){
    DEBUG("Creating configuration.");
    Manager::_configuration = new ConfigurationExternal(_p);
    DEBUG("Configuration created.");
    // For instrumented application we do not care if synchronous of not
    // (we count iterations).
    _p.synchronousWorkers = false;
}

ManagerInstrumented::ManagerInstrumented(nn::socket& riffSocket,
                                 int chid,
                                 Parameters nornirParameters):
            Manager(nornirParameters), _monitor(riffSocket, chid){
    Manager::_configuration = new ConfigurationExternal(_p);
    // For instrumented application we do not care if synchronous of not (we
    // count iterations).
    _p.synchronousWorkers = false;
}

ManagerInstrumented::~ManagerInstrumented(){
    if(Manager::_configuration){
        delete Manager::_configuration;
    }
    if(_selector){
        delete _selector;
    }
}

void ManagerInstrumented::waitForStart(){
    Manager::_pid = _monitor.waitStart();
    dynamic_cast<KnobVirtualCores*>(_configuration->getKnob(KNOB_VIRTUAL_CORES))->changeMax(_monitor.getTotalThreads());
    dynamic_cast<KnobMappingExternal*>(_configuration->getKnob(KNOB_MAPPING))->setPid(_pid);
    dynamic_cast<KnobClkModEmulated*>(_configuration->getKnob(KNOB_CLKMOD_EMULATED))->setPid(_pid);
}

MonitoredSample ManagerInstrumented::getSample(){
    MonitoredSample sample;
    if(!_monitor.getSample(sample)){
        _terminated = true;
    }
    // Knarr may return inconsistent data for latency and
    // utilization factor when performs sampling.
    // Check if this is the case.
    if(sample.inconsistent){
        if(_p.requirements.minUtilization != NORNIR_REQUIREMENT_UNDEF ||
           _p.requirements.maxUtilization != NORNIR_REQUIREMENT_UNDEF){
            throw std::runtime_error("You specified requirements on loadPercentage but instrumenter is "
                                     "providing inconsistent loadPercentage values. Please call "
                                     "setConfiguration with samplingLengthMs = 0 on nornir's Instrumenter "
                                     "to fix this issue.");
        }
    }
    if(sample.inconsistent){
        if(_p.requirements.latency != NORNIR_REQUIREMENT_UNDEF){
            throw std::runtime_error("You specified requirements on latency but instrumenter is "
                                     "providing inconsistent latency values. Please call "
                                     "setConfiguration with samplingLengthMs = 0 on nornir's Instrumenter "
                                     "to fix this issue.");
        }
    }
    return sample;
}

MonitoredSample ManagerInstrumented::clearStoredSample(){
    return getSample();
}

ulong ManagerInstrumented::getExecutionTime(){
    return _monitor.getExecutionTime();
}

void ManagerInstrumented::terminationManagement(){
    // We need to overwrite the tasks count because
    // due to sampling, the last batch of tasks may still
    // not have been communicated to the manager.
    // By doing so, we are sure that _totalTasks
    // represents the total amount of processed tasks.
    _totalTasks = _monitor.getTotalTasks();
}

void ManagerInstrumented::shrinkPause(){
    kill(_pid, SIGSTOP);
}

void ManagerInstrumented::stretchPause(){
    kill(_pid, SIGCONT);
}

ManagerBlackBox::ManagerBlackBox(pid_t pid, Parameters nornirParameters):
        Manager(nornirParameters), _process(nornirParameters.mammut.getInstanceTask()->getProcessHandler(pid)){
    Manager::_pid = pid;
    Manager::_configuration = new ConfigurationExternal(_p);
    // For blackbox application we do not care if synchronous of not
    // (we count instructions).
    _p.synchronousWorkers = false;
    // Check supported requirements.
    if(isPrimaryRequirement(_p.requirements.throughput) ||
       _p.requirements.maxUtilization != NORNIR_REQUIREMENT_UNDEF ||
       isPrimaryRequirement(_p.requirements.executionTime) ||
       _p.requirements.latency != NORNIR_REQUIREMENT_UNDEF){
        throw std::runtime_error("ManagerBlackBox. Unsupported requirement.");
    }
    _startTime = 0;
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
    // However if we have been required to monitor only the ROI, we
    // wait for ROI start by waiting for the roiFile creation.
    if(_p.roiFile.compare("")){
        while(!mammut::utils::existsFile(_p.roiFile)){
            usleep(1000);
        }
    }
    _startTime = getMillisecondsTime();
    dynamic_cast<KnobMappingExternal*>(_configuration->getKnob(KNOB_MAPPING))->setProcessHandler(_process);
    dynamic_cast<KnobClkModEmulated*>(_configuration->getKnob(KNOB_CLKMOD_EMULATED))->setProcessHandler(_process);
    _process->resetInstructions(); // To remove those executed before entering ROI
}

static bool isRunning(pid_t pid) {
    return !kill(pid, 0);
}

MonitoredSample ManagerBlackBox::getSample(){
    MonitoredSample sample;
    double instructions = 0;
    if(!isRunning(_process->getId()) ||
       (_p.roiFile.compare("") && !mammut::utils::existsFile(_p.roiFile))){
        _terminated = true;
        return sample;
    }
    _process->getAndResetInstructions(instructions);
    sample.throughput = instructions / ((getMillisecondsTime() - _lastStoredSampleMs) / 1000.0);
    sample.latency = -1; // Not used.
    sample.loadPercentage = 100.0; // We do not know what's the input bandwidth.
    sample.numTasks = instructions; // We consider a task to be an instruction.
    return sample;
}

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

void ManagerFastFlow::waitForStart(){
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

void ManagerFastFlow::askForSample(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->askForSample();
    }
}

MonitoredSample ManagerFastFlow::getSampleResponse(){
    MonitoredSample sample;
    uint numActiveWorkers = _activeWorkers.size();
    for(size_t i = 0; i < numActiveWorkers; i++){
        MonitoredSample tmp;
        AdaptiveNode* w = _activeWorkers.at(i);
        w->getSampleResponse(tmp, _samples->average().latency);
        sample.loadPercentage += tmp.loadPercentage;
        sample.numTasks += tmp.numTasks;
        sample.latency += tmp.latency;
        sample.throughput += tmp.throughput;
    }
    sample.loadPercentage /= numActiveWorkers;
    sample.latency /= numActiveWorkers;
    return sample;
}

MonitoredSample ManagerFastFlow::getSample(){
    askForSample();
    return getSampleResponse();
}

ManagerFastFlow::ManagerFastFlow(ff_farm<>* farm,
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
}

ManagerFastFlow::~ManagerFastFlow(){
    delete _samples;
    delete _variations;
    if(_selector){
        delete _selector;
    }
    if(Manager::_configuration){
        delete Manager::_configuration;
    }
}

void ManagerFastFlow::initNodesPreRun() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->initPreRun(_p, NODE_TYPE_WORKER, &_terminated);
    }
    if (_emitter) {
        _emitter->initPreRun(_p, NODE_TYPE_EMITTER, &_terminated, _farm->getlb());
    } else {
        throw runtime_error("Emitter is needed to use the manager.");
    }
    if (_collector) {
        _collector->initPreRun(_p, NODE_TYPE_COLLECTOR, &_terminated,
                               _farm->getgt());
    }
}

void ManagerFastFlow::initNodesPostRun() {
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

void ManagerFastFlow::postConfigurationManagement(){
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

void ManagerFastFlow::terminationManagement(){
    DEBUG("Terminating...wait freezing.");
    _farm->wait_freezing();
    _farm->wait();
    DEBUG("Terminated.");
}

ulong ManagerFastFlow::getExecutionTime(){
    return _farm->ffTime();
}

void ManagerFastFlow::shrinkPause(){
    KnobVirtualCoresFarm* k = dynamic_cast<KnobVirtualCoresFarm*>(_configuration->getKnob(KNOB_VIRTUAL_CORES));
    k->prepareToFreeze();
    k->freeze();
}

void ManagerFastFlow::stretchPause(){
    KnobVirtualCoresFarm* k = dynamic_cast<KnobVirtualCoresFarm*>(_configuration->getKnob(KNOB_VIRTUAL_CORES));
    size_t v = k->getRealValue();
    k->prepareToRun(v);
    k->run(v);
}

}
