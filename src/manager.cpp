/*
 * manager.cpp
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

#include "./manager.hpp"
#include "parameters.hpp"
#include "predictors.hpp"
#include "./node.hpp"
#include "utils.hpp"

#include <ff/farm.hpp>
#include <mammut/module.hpp>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

#include <cmath>
#include <iostream>
#include <limits>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_MANAGER
#define DEBUG(x) do { cerr << "[Manager] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace adpff{

class Parameters;
class ManagerFarm;

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
double ManagerFarm<lb_t, gt_t>::getMaxPredictionErrorPrimary() const{
    double r = _p.maxPrimaryPredictionError;
    if(_p.strategyPredictionErrorPrimary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getPrimaryValue(_samples->coefficientVariation()));
    }
    return r;
}

template <typename lb_t, typename gt_t>
double ManagerFarm<lb_t, gt_t>::getMaxPredictionErrorSecondary() const{
    double r = _p.maxSecondaryPredictionError;
    if(_p.strategyPredictionErrorSecondary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getSecondaryValue(_samples->coefficientVariation()));
    }
    return r;
}

template <typename lb_t, typename gt_t>
double ManagerFarm<lb_t, gt_t>::getPrimaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            return sample.utilization;
        }break;
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
double ManagerFarm<lb_t, gt_t>::getPrimaryValue() const{
    return getPrimaryValue(_samples->average());
}

template <typename lb_t, typename gt_t>
double ManagerFarm<lb_t, gt_t>::getSecondaryValue() const{
    return getSecondaryValue(_samples->average());
}

template <typename lb_t, typename gt_t>
bool ManagerFarm<lb_t, gt_t>::isContractViolated() const{
    double tolerance = 0;
    double maxError = getMaxPredictionErrorPrimary();
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            tolerance = ((_p.overloadThresholdFarm -
                          _p.underloadThresholdFarm) * maxError) / 100.0;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            tolerance = (_p.requiredBandwidth * maxError) / 100.0;
        }break;
        case CONTRACT_POWER_BUDGET:{
            tolerance = (_p.powerBudget * maxError) / 100.0;
        }break;
        default:{
            return false;
        }break;
    }
    return !isFeasiblePrimaryValue(getPrimaryValue(), tolerance);
}

template <typename lb_t, typename gt_t>
bool ManagerFarm<lb_t, gt_t>::terminated(){
    /**
     * We do not need to wait if the emitter is terminated.
     * Indeed, if the workers terminated, the emitter surely terminated
     * too.
     */

    for(size_t i = 0; i < _activeWorkers.size(); i++){
        if(!_activeWorkers.at(i)->isTerminated()){
            return false;
        }else{
            DEBUG("Worker " << i << " terminated.");
        }
    }

    if(_collector &&
       !_collector->isTerminated()){
        return false;
    }else{
        DEBUG("Collector terminated.");
    }

    return true;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::changeKnobs(){
    KnobsValues values = _calibrator->getNextKnobsValues();
    if(values != _configuration.getRealValues()){
        _configuration.setValues(values);
        _activeWorkers = dynamic_cast<const KnobWorkers*>(_configuration.getKnob(KNOB_TYPE_WORKERS))->getActiveWorkers();

        /****************** Clean state ******************/
        _lastStoredSampleMs = getMillisecondsTime();
        _samples->reset();
        _counter->reset();
        _totalTasks = 0;
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::observe(){
    if(_p.observer){
        const KnobMapping* kMapping = dynamic_cast<const KnobMapping*>(_configuration.getKnob(KNOB_TYPE_MAPPING));
        _p.observer->observe(_lastStoredSampleMs,
                             _configuration.getRealValue(KNOB_TYPE_WORKERS),
                             _configuration.getRealValue(KNOB_TYPE_FREQUENCY),
                             kMapping->getEmitterVirtualCore(),
                             kMapping->getWorkersVirtualCore(),
                             kMapping->getCollectorVirtualCore(),
                             _samples->getLastSample().bandwidth,
                             _samples->average().bandwidth,
                             _samples->coefficientVariation().bandwidth,
                             _samples->average().utilization,
                             _samples->average().watts);
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
        w->getSampleResponse(tmp, _p.strategyPolling,
                             _samples->average().latency);
        sample += tmp;
    }
    sample.loadPercentage /= numActiveWorkers;
    sample.latency /= numActiveWorkers;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::storeNewSample(){
    MonitoredSample sample;
    WorkerSample ws;
    Joules joules;

    askForWorkersSamples();
    getWorkersSamples(ws);

    _totalTasks += ws.tasksCount;
    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        if(_remainingTasks > ws.tasksCount){
            _remainingTasks -= ws.tasksCount;
        }else{
            _remainingTasks = 0;
        }
    }

    switch(_counter->getType()){
        case COUNTER_CPUS:{
            joules = ((CounterCpus*) _counter)->getJoulesCoresAll();
        }break;
        default:{
            joules = _counter->getJoules();
        }break;
    }

    double now = getMillisecondsTime();
    double durationSecs = (now - _lastStoredSampleMs) / 1000.0;
    _lastStoredSampleMs = now;

    sample.watts = joules / durationSecs;
    sample.utilization = ws.loadPercentage;
    // ATTENTION: Bandwidth is not the number of task since the
    //            last observation but the number of expected
    //            tasks that will be processed in 1 second.
    //            For this reason, if we sum all the bandwidths in
    //            the result observation file, we may have an higher
    //            number than the number of tasks.
    sample.bandwidth = ws.bandwidthTotal;
    sample.latency = ws.latency;

    _counter->reset();
    _samples->add(sample);

    DEBUGB(samplesFile << *_samples << "\n");
}

template <typename lb_t, typename gt_t>
bool ManagerFarm<lb_t, gt_t>::persist() const{
    bool r = false;
    switch(_p.strategyPersistence){
        case STRATEGY_PERSISTENCE_SAMPLES:{
            r = _samples->size() < _p.persistenceValue;
        }break;
        case STRATEGY_PERSISTENCE_TASKS:{
            r = _totalTasks < _p.persistenceValue;
        }break;
        case STRATEGY_PERSISTENCE_VARIATION:{
            const MonitoredSample& variation =
                    _samples->coefficientVariation();
            r = getPrimaryValue(variation) < _p.persistenceValue &&
                getSecondaryValue(variation) < _p.persistenceValue;
        }break;
    }
    return r;
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initPredictors(){
    if(_p.strategyCalibration == STRATEGY_CALIBRATION_RANDOM){
        ; //CREARE
    }else{
        _calibrator = new CalibratorLowDiscrepancy(p, *this);
    }
}

#include <string>
Parameters& validate(Parameters& p){
    ParametersValidation apv = p.validate();
    if(apv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + std::to_string(apv));
    }
    return p;
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::ManagerFarm(ff_farm<lb_t, gt_t>* farm, Parameters parameters):
        _farm(farm),
        _p(validate(parameters)),
        _startTimeMs(0),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _counter(_p.mammut.getInstanceEnergy()->getCounter()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _numCpus(_topology->getCpus().size()),
        _numPhysicalCores(_topology->getPhysicalCores().size()),
        _numPhysicalCoresPerCpu(_topology->getCpu(0)->getPhysicalCores().
                                size()),
        _numVirtualCoresPerPhysicalCore(_topology->getPhysicalCore(0)->
                                        getVirtualCores().size()),
        _emitter(dynamic_cast<AdaptiveNode*>(farm->getEmitter())),
        _collector(dynamic_cast<AdaptiveNode*>(farm->getCollector())),
        _activeWorkers(convertWorkers(farm->getWorkers())),
        _configuration<lb_t, gt_t>(_p, *farm),
        _samples(NULL),
        _totalTasks(0),
        _remainingTasks(0),
        _deadline(0),
        _lastStoredSampleMs(0),
        _calibrator(NULL),
        _primaryPredictor(NULL),
        _secondaryPredictor(NULL),
        _primaryPrediction(0),
        _secondaryPrediction(0){

    _samples = NULL;
    switch(_p.strategySmoothing){
    case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
        _samples = new MovingAverageSimple<MonitoredSample>(_p.smoothingFactor);
    }break;
    case STRATEGY_SMOOTHING_EXPONENTIAL:{
        _samples = new MovingAverageExponential<MonitoredSample>(_p.smoothingFactor);
    }break;
    }

    DEBUGB(samplesFile.open("samples.csv"));
}

template <typename lb_t, typename gt_t>
ManagerFarm<lb_t, gt_t>::~ManagerFarm(){
    delete _samples;
    if(_primaryPredictor){
        delete _primaryPredictor;
    }
    if(_secondaryPredictor){
        delete _secondaryPredictor;
    }
    if(_calibrator){
        delete _calibrator;
    }
    DEBUGB(samplesFile.close());
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initNodesPreRun() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->initPreRun(_p.mammut, _p.archData.ticksPerNs, NODE_TYPE_WORKER);
    }
    if (_emitter) {
        _emitter->initPreRun(_p.mammut, _p.archData.ticksPerNs, NODE_TYPE_EMITTER);
    } else {
        throw runtime_error("Emitter is needed to use the manager.");
    }
    if (_collector) {
        _collector->initPreRun(_p.mammut, _p.archData.ticksPerNs, NODE_TYPE_COLLECTOR);
    }
}

template <typename lb_t, typename gt_t>
void ManagerFarm<lb_t, gt_t>::initNodesPostRun() {
    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->initPostRun();
    }
    _emitter->initPostRun();
    if (_collector) {
        _collector->initPostRun();
    }
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
void ManagerFarm<lb_t, gt_t>::run(){
    initNodesPreRun();
    _farm->run_then_freeze();
    initNodesPostRun();
    _configuration.maxAllKnobs();

    _startTimeMs = getMillisecondsTime();
    _counter->reset();
    _lastStoredSampleMs = _startTimeMs;
    if(_p.observer){
        _p.observer->_startMonitoringMs = _lastStoredSampleMs;
    }

    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
        _remainingTasks = _p.expectedTasksNumber;
        _deadline = getMillisecondsTime()/1000.0 +
                    _p.requiredCompletionTime;
    }

    initPredictors();

    double microsecsSleep = 0;
    if(_p.contractType == CONTRACT_NONE){
        _farm->wait();
        storeNewSample();
        observe();
    }else{
        /* Force the first calibration point. **/
        assert(_calibrator);
        changeKnobs();

        double startSample = getMillisecondsTime();

        while(!terminated()){
            double overheadMs = getMillisecondsTime() - startSample;
            microsecsSleep = ((double)_p.samplingInterval - overheadMs)*
                              (double)MAMMUT_MICROSECS_IN_MILLISEC;
            if(microsecsSleep < 0){
                microsecsSleep = 0;
            }
            usleep(microsecsSleep);
            startSample = getMillisecondsTime();

            storeNewSample();
            DEBUG("New sample stored.");

            if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                uint now = getMillisecondsTime()/1000.0;
                if(now >= _deadline){
                    _p.requiredBandwidth = numeric_limits<double>::max();
                }else{
                    _p.requiredBandwidth = _remainingTasks /
                                           (_deadline - now);
                }
            }

            observe();

            if(!persist()){
                assert(_calibrator);
                changeKnobs();
                startSample = getMillisecondsTime();
            }
        }
        DEBUG("Terminated.");
    }

    uint duration = getMillisecondsTime() - _startTimeMs;
    if(_p.observer){
        vector<CalibrationStats> cs;
        if(_calibrator){
            cs = _calibrator->getCalibrationsStats();
            _p.observer->calibrationStats(cs, duration);
        }
        _p.observer->summaryStats(cs, duration);
    }

    cleanNodes();
}

double Observer::calibrationDurationToPerc(const CalibrationStats& cs,
                                 uint durationMs){
    return ((double)cs.duration /
            (double)durationMs) * 100.0;
}

Observer::Observer(string statsFile, string calibrationFile, string summaryFile):
        _startMonitoringMs(0),
        _totalWatts(0),
        _totalBw(0),
        _numSamples(0){
    _statsFile.open(statsFile.c_str());
    _calibrationFile.open(calibrationFile.c_str());
    _summaryFile.open(summaryFile.c_str());
    if(!_statsFile.is_open() ||
       !_calibrationFile.is_open() ||
       !_summaryFile.is_open()){
        throw runtime_error("Observer: Impossible to open file.");
    }
    _statsFile << "TimestampMillisecs" << "\t";
    _statsFile << "[[EmitterVc][WorkersVc][CollectorVc]]" << "\t";
    _statsFile << "Workers" << "\t";
    _statsFile << "Frequency" << "\t";
    _statsFile << "CurrentBandwidth" << "\t";
    _statsFile << "SmoothedBandwidth" << "\t";
    _statsFile << "CoeffVarBandwidth" << "\t";
    _statsFile << "SmoothedUtilization" << "\t";
    _statsFile << "SmoothedWatts" << "\t";
    _statsFile << endl;

    _calibrationFile << "NumSteps" << "\t";
    _calibrationFile << "Duration" << "\t";
    _calibrationFile << "Time%" << "\t";
    _calibrationFile << endl;

    _summaryFile << "Watts" << "\t";
    _summaryFile << "Bandwidth" << "\t";
    _summaryFile << "CompletionTime" << "\t";
    _summaryFile << "Calibration%" << "\t";
    _summaryFile << endl;
}

Observer::~Observer(){
    _statsFile.close();
    _calibrationFile.close();
    _summaryFile.close();
}

void Observer::observe(unsigned int timeStamp,
                     size_t workers,
                     Frequency frequency,
                     const VirtualCore* emitterVirtualCore,
                     const vector<VirtualCore*>& workersVirtualCore,
                     const VirtualCore* collectorVirtualCore,
                     double currentBandwidth,
                     double smoothedBandwidth,
                     double coeffVarBandwidth,
                     double smoothedUtilization,
                     Joules smoothedWatts){
    _statsFile << timeStamp - _startMonitoringMs << "\t";
    _statsFile << "[";
    if(emitterVirtualCore){
        _statsFile << "[" << emitterVirtualCore->getVirtualCoreId() << "]";
    }

    _statsFile << "[";
    for(size_t i = 0; i < workersVirtualCore.size(); i++){
        _statsFile << workersVirtualCore.at(i)->getVirtualCoreId() << ",";
    }
    _statsFile << "]";

    if(collectorVirtualCore){
        _statsFile << "[" << collectorVirtualCore->getVirtualCoreId() << "]";
    }
    _statsFile << "]" << "\t";

    _statsFile << workers << "\t";
    _statsFile << frequency << "\t";
    _statsFile << currentBandwidth << "\t";
    _statsFile << smoothedBandwidth << "\t";
    _statsFile << coeffVarBandwidth << "\t";
    _statsFile << smoothedUtilization << "\t";

    _statsFile << smoothedWatts << "\t";

    _statsFile << endl;

    _totalWatts += smoothedWatts;
    _totalBw += currentBandwidth;
    _numSamples++;
}

void Observer::calibrationStats(const vector<CalibrationStats>&
                              calibrationStats,
                              uint durationMs){

    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        _calibrationFile << cs.numSteps << "\t";
        _calibrationFile << cs.duration << "\t";
        _calibrationFile << calibrationDurationToPerc(cs, durationMs) << "\t";
        _calibrationFile << endl;
    }
}

void Observer::summaryStats(const vector<CalibrationStats>&
                          calibrationStats,
                          uint durationMs){
    double totalCalibrationPerc = 0.0;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        totalCalibrationPerc += calibrationDurationToPerc(cs, durationMs);
    }

    _summaryFile << _totalWatts / (double) _numSamples << "\t";
    _summaryFile << _totalBw / (double) _numSamples << "\t";
    _summaryFile << (double) durationMs / 1000.0 << "\t";
    _summaryFile << totalCalibrationPerc << "\t";
    _summaryFile << endl;
}

}
