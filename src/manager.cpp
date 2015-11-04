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

void ManagerFarm::setDomainToHighestFrequency(const Domain* domain){
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

double ManagerFarm::getMaxPredictionErrorPrimary() const{
    double r = _p.maxPrimaryPredictionError;
    if(_p.strategyPredictionErrorPrimary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getPrimaryValue(_samples->coefficientVariation()));
    }
    return r;
}

double ManagerFarm::getMaxPredictionErrorSecondary() const{
    double r = _p.maxSecondaryPredictionError;
    if(_p.strategyPredictionErrorSecondary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getSecondaryValue(_samples->coefficientVariation()));
    }
    return r;
}

double ManagerFarm::getPrimaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            return sample.utilization;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return sample.bandwidth;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return sample.watts.cores;
        }break;
        default:{
            return 0;
        }break;
    }
}

double ManagerFarm::getSecondaryValue(const MonitoredSample& sample) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return sample.watts.cores;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return sample.bandwidth;
        }break;
        default:{
            return 0;
        }break;
    }
}

double ManagerFarm::getPrimaryValue() const{
    return getPrimaryValue(_samples->average());
}

double ManagerFarm::getSecondaryValue() const{
    return getSecondaryValue(_samples->average());
}

bool ManagerFarm::isContractViolated() const{
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

bool ManagerFarm::isFeasiblePrimaryValue(double value, double tolerance) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            return value > _p.underloadThresholdFarm - tolerance &&
                   value < _p.overloadThresholdFarm + tolerance;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            return value > _p.requiredBandwidth - tolerance;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return value < _p.powerBudget + tolerance;
        }break;
        default:{
            return false;
        }break;
    }
    return false;
}

//TODO Move in predictors.cpp
double ManagerFarm::getVoltage(const KnobsValues& values) const{
    VoltageTableKey key(values[KNOB_TYPE_WORKERS], values[KNOB_TYPE_FREQUENCY]);
    VoltageTableIterator it = _voltageTable.find(key);
    if(it != _voltageTable.end()){
        return it->second;
    }else{
        throw runtime_error("Frequency and/or number of virtual cores "
                                 "not found in voltage table.");
    }
}

bool ManagerFarm::isBestSuboptimalValue(double x, double y) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            // Concerning utilization factors, if both are suboptimal,
            // we prefer the closest to the lower bound.
            double distanceX, distanceY;
            distanceX = _p.underloadThresholdFarm - x;
            distanceY = _p.underloadThresholdFarm - y;
            if(distanceX > 0 && distanceY < 0){
                return true;
            }else if(distanceX < 0 && distanceY > 0){
                return false;
            }else{
                return abs(distanceX) < abs(distanceY);
            }
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            // Concerning bandwidths, if both are suboptimal,
            // we prefer the higher one.
            return x > y;
        }break;
        case CONTRACT_POWER_BUDGET:{
            // Concerning power budgets, if both are suboptimal,
            // we prefer the lowest one.
            return x < y;
        }break;
        default:{
            ;
        }break;
    }
    return false;
}

bool ManagerFarm::isBestSecondaryValue(double x, double y) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_BANDWIDTH:{
            return x < y;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return x > y;
        }break;
        default:{
            ;
        }break;
    }
    return false;
}

KnobsValues ManagerFarm::getNewKnobsValues(){
    KnobsValues bestValues;
    KnobsValues bestSuboptimalValues = _configuration.getRealValues();

    double primaryPrediction = 0;
    double secondaryPrediction = 0;

    double bestPrimaryPrediction = 0;
    double bestSecondaryPrediction = 0;
    double bestSuboptimalValue = getPrimaryValue();

    bool feasibleSolutionFound = false;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            // We have to minimize the power/energy.
            bestSecondaryPrediction = numeric_limits<double>::max();
        }break;
        case CONTRACT_POWER_BUDGET:{
            // We have to maximize the bandwidth.
            bestSecondaryPrediction = numeric_limits<double>::min();
        }break;
        default:{
            ;
        }break;
    }

    _primaryPredictor->prepareForPredictions();
    _secondaryPredictor->prepareForPredictions();

    unsigned int remainingTime = 0;
    vector<KnobsValues> combinations = _configuration.getAllRealCombinations();
    for(size_t i = 0; i < combinations.size(); i++){
        KnobsValues currentValues = combinations.at(i);
        primaryPrediction = _primaryPredictor->predict(currentValues);
        switch(_p.contractType){
            case CONTRACT_PERF_COMPLETION_TIME:{
                remainingTime = (double) _remainingTasks /
                                 primaryPrediction;
            }break;
            case CONTRACT_PERF_UTILIZATION:{
                primaryPrediction = (_samples->average().bandwidth /
                                     primaryPrediction) *
                                     _samples->average().utilization;
            }break;
            default:{
                ;
            }
        }

        if(isFeasiblePrimaryValue(primaryPrediction)){
            secondaryPrediction = _secondaryPredictor->
                                  predict(currentValues);
            if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                secondaryPrediction *= remainingTime;
            }
            if(isBestSecondaryValue(secondaryPrediction,
                                    bestSecondaryPrediction)){
                bestValues = currentValues;
                feasibleSolutionFound = true;
                bestPrimaryPrediction = primaryPrediction;
                bestSecondaryPrediction = secondaryPrediction;
            }
        }else if(!feasibleSolutionFound &&
                 isBestSuboptimalValue(primaryPrediction,
                                       bestSuboptimalValue)){
            bestSuboptimalValue = primaryPrediction;
            bestSuboptimalValues = currentValues;
        }
    }

    if(feasibleSolutionFound){
        _primaryPrediction = bestPrimaryPrediction;
        _secondaryPrediction = bestSecondaryPrediction;
        return bestValues;
    }else{
        _primaryPrediction = bestSuboptimalValue;
        // TODO: This check now works because both service time and power are always  > 0
        // In the future we must find another way to indicat that secondary prediction
        // has not been done.
        _secondaryPrediction = -1;
        return bestSuboptimalValues;
    }
}

bool ManagerFarm::terminated(){
    if(_emitter &&
       _emitter->isTerminated()){
        return true;
    }

    for(size_t i = 0; i < _activeWorkers.size(); i++){
        if(_activeWorkers.at(i)->isTerminated()){
            return true;
        }
    }

    if(_collector &&
       _collector->isTerminated()){
        return true;
    }

    return false;
}

void ManagerFarm::changeRelative(KnobsValues values){
    _configuration.setRelativeValues(values);
    _activeWorkers = dynamic_cast<const KnobWorkers*>(_configuration.getKnob(KNOB_TYPE_WORKERS))->getActiveWorkers();

    /****************** Clean state ******************/
    _lastStoredSampleMs = getMillisecondsTime();
    _samples->reset();
    _energy->resetCountersCpu();
    _totalTasks = 0;
}

void ManagerFarm::observe(){
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

void ManagerFarm::askForWorkersSamples(){
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->askForSample();
    }
}

void ManagerFarm::getWorkersSamples(WorkerSample& sample){
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

void ManagerFarm::storeNewSample(){
    MonitoredSample sample;
    WorkerSample ws;
    JoulesCpu joules;

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

    vector<CounterCpu*> energyCounters = _energy->getCountersCpu();
    for(size_t i = 0; i < energyCounters.size(); i++){
        joules += energyCounters.at(i)->getJoules();
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

    _energy->resetCountersCpu();
    _samples->add(sample);

    DEBUGB(samplesFile << *_samples << "\n");
}

bool ManagerFarm::persist() const{
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

void ManagerFarm::initCalibrator(){
    if(_p.strategyCalibration == STRATEGY_CALIBRATION_RANDOM){
        ;
    }else{
        _calibrator = new CalibratorLowDiscrepancy(*this);
    }
}

void ManagerFarm::initPredictors(){
    PredictorType primary, secondary;
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            primary = PREDICTION_BANDWIDTH;
            secondary = PREDICTION_POWER;
        }break;
        case CONTRACT_POWER_BUDGET:{
            primary = PREDICTION_POWER;
            secondary = PREDICTION_BANDWIDTH;
        }break;
        default:{
            return;
        }break;
    }

    switch(_p.strategyPrediction){
        case STRATEGY_PREDICTION_SIMPLE:{
            _primaryPredictor = new PredictorSimple(primary, *this);
            _secondaryPredictor = new PredictorSimple(secondary, *this);
            _calibrator = NULL;
        }break;
        case STRATEGY_PREDICTION_REGRESSION_LINEAR:{
            _primaryPredictor = new PredictorLinearRegression(primary,
                                                              *this);
            _secondaryPredictor = new PredictorLinearRegression(secondary,
                                                                *this);
            initCalibrator();
        }break;
        default:{
            ;
        }break;
    }
}

Parameters& validate(Parameters& p){
    ParametersValidation apv = p.validate();
    if(apv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + apv);
    }
    return p;
}

ManagerFarm::ManagerFarm(ff_farm<>* farm, Parameters parameters):
        _farm(farm),
        _p(validate(parameters)),
        _startTimeMs(0),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _energy(_p.mammut.getInstanceEnergy()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _numCpus(_topology->getCpus().size()),
        _numPhysicalCores(_topology->getPhysicalCores().size()),
        _numPhysicalCoresPerCpu(_topology->getCpu(0)->getPhysicalCores().
                                size()),
        _numVirtualCoresPerPhysicalCore(_topology->getPhysicalCore(0)->
                                        getVirtualCores().size()),
        _configuration(_p, *farm){

    /** If voltage table file is specified, then load the table. **/
    if(_p.archData.voltageTableFile.compare("")){
        loadVoltageTable(_voltageTable,
                         _p.archData.voltageTableFile);
    }

    _samples = NULL;
    switch(_p.strategySmoothing){
    case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
        _samples = new MovingAverageSimple<MonitoredSample>
                                (_p.smoothingFactor);
    }break;
    case STRATEGY_SMOOTHING_EXPONENTIAL:{
        _samples = new MovingAverageExponential<MonitoredSample>
                                (_p.smoothingFactor);
    }break;
    }

    DEBUGB(samplesFile.open("samples.csv"));
}

ManagerFarm::~ManagerFarm(){
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

void ManagerFarm::initNodes() {
    _emitter = dynamic_cast<AdaptiveNode*>(_farm->getEmitter());
    _collector = dynamic_cast<AdaptiveNode*>(_farm->getCollector());
    svector<ff_node*> w = _farm->getWorkers();
    for(size_t i = 0; i < w.size(); i++){
        _activeWorkers.push_back(dynamic_cast<AdaptiveNode*>(w[i]));
    }

    for (size_t i = 0; i < _activeWorkers.size(); i++) {
        _activeWorkers.at(i)->init(_p.mammut, _p.archData.ticksPerNs);
    }
    if (_emitter) {
        _emitter->init(_p.mammut, _p.archData.ticksPerNs);
    } else {
        throw runtime_error("Emitter is needed to use the manager.");
    }
    if (_collector) {
        _collector->init(_p.mammut, _p.archData.ticksPerNs);
    }
}

void ManagerFarm::cleanNodes() {
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

void ManagerFarm::run(){
    _farm->run_then_freeze(_farm->getNWorkers());

    initNodes();
    _configuration.maxAllKnobs();

    _startTimeMs = getMillisecondsTime();
    _energy->resetCountersCpu();
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
        if(_calibrator){
            changeRelative(_calibrator->getNextKnobsValues());
        }

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
                bool reconfigurationRequired = false;
                KnobsValues nextValues;

                if(_calibrator){
                    nextValues = _calibrator->getNextKnobsValues();
                    reconfigurationRequired = true;
                }else if(isContractViolated()){
                    nextValues = getNewKnobsValues();
                    reconfigurationRequired = true;
                }

                if(reconfigurationRequired){
                    changeRelative(nextValues);
                    startSample = getMillisecondsTime();
                }
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
    _statsFile << "SmoothedWattsCpu" << "\t";
    _statsFile << "SmoothedWattsCores" << "\t";
    _statsFile << "SmoothedWattsGraphic" << "\t";
    _statsFile << "SmoothedWattsDram" << "\t";
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
                     JoulesCpu smoothedWatts){
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

    _statsFile << smoothedWatts.cpu << "\t";
    _statsFile << smoothedWatts.cores << "\t";
    _statsFile << smoothedWatts.graphic << "\t";
    _statsFile << smoothedWatts.dram << "\t";

    _statsFile << endl;

    _totalWatts += smoothedWatts.cores;
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

FarmConfiguration::FarmConfiguration(const Parameters& p, ff::ff_farm<>& farm):_p(p){
    _knobs[KNOB_TYPE_WORKERS] = new KnobWorkers(p.knobWorkers, farm);
    _knobs[KNOB_TYPE_MAPPING] = new KnobMapping(p.knobMapping,
                                                p.knobMappingEmitter,
                                                p.knobMappingCollector,
                                                p.knobHyperthreading,
                                                p.mammut,
                                                dynamic_cast<AdaptiveNode*>(farm.getEmitter()),
                                                dynamic_cast<AdaptiveNode*>(farm.getCollector()),
                                                *((KnobWorkers*)_knobs[KNOB_TYPE_WORKERS]));
    _knobs[KNOB_TYPE_FREQUENCY] = new KnobFrequency(p.knobFrequencies,
                                                    p.mammut,
                                                    p.turboBoost,
                                                    p.strategyInactiveVirtualCores,
                                                    p.strategyUnusedVirtualCores,
                                                    *((KnobMapping*)_knobs[KNOB_TYPE_MAPPING]));

    std::vector<std::vector<double>> values;
    std::vector<double> accum;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        values.push_back(_knobs[i]->getAllowedValues());
    }
    combinations(values, 0, accum);
}

FarmConfiguration::~FarmConfiguration(){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        delete _knobs[i];
    }
}

//TODO: Works even if a vector is empty? (i.e. a knob has no values)
void FarmConfiguration::combinations(vector<vector<double> > array, size_t i, vector<double> accum){
    if(i == array.size()){
        KnobsValues kv;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            kv[(KnobType) i] = accum.at(i);
        }
        _combinations.push_back(kv);
    }else{
        vector<double> row = array.at(i);
        for(size_t j = 0; j < row.size(); ++j){
            vector<double> tmp(accum);
            tmp.push_back(row[j]);
            combinations(array, i+1, tmp);
        }
    }
}

const std::vector<KnobsValues>& FarmConfiguration::getAllRealCombinations(){
    return _combinations;
}

void FarmConfiguration::setFastReconfiguration(){
    if(_p.fastReconfiguration){
        ((KnobFrequency*) _knobs[KNOB_TYPE_FREQUENCY])->setRelativeValue(100.0);
    }
}

const Knob* FarmConfiguration::getKnob(KnobType t) const{
    return _knobs[t];
}

void FarmConfiguration::maxAllKnobs(){
    DEBUG("Maxing all the knobs.");
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[(KnobType) i]->setToMax();
    }
}

double FarmConfiguration::getRealValue(KnobType t) const{
    return _knobs[t]->getRealValue();
}

KnobsValues FarmConfiguration::getRealValues() const{
    KnobsValues kv;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        kv[(KnobType) i] = getRealValue((KnobType) i);
    }
    return kv;
}

double FarmConfiguration::getRelativeValue(KnobType t) const{
    return _knobs[t]->getRelativeValue();
}

void FarmConfiguration::setRelativeValues(const KnobsValues& values){
    // Fast reconfiguration is valid only for knobs changed before
    // the frequency knob.
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[i]->setRelativeValue(values[(KnobType)i]);
    }
    DEBUG("Changed relative knobs values.");
}

void FarmConfiguration::setRealValues(const KnobsValues& values){
    // Fast reconfiguration is valid only for knobs changed before
    // the frequency knob.
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[i]->setRealValue(values[(KnobType)i]);
    }
    DEBUG("Changed real knobs values.");
}

}
