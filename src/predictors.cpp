/*
 * predictors.cpp
 *
 * Created on: 10/07/2015
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

/*!
 * \file predictors.cpp
 * \brief Predictors used by adaptive farm.
 */

#include "predictors.hpp"
#include "manager.hpp"

#include <mammut/cpufreq/cpufreq.hpp>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_PREDICTORS
#define DEBUG(x) do { std::cerr << "[Predictors] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

#define NO_PREDICTION DBL_MIN

namespace adpff{

using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

static double getVoltage(VoltageTable table, uint workers, Frequency frequency){
    VoltageTableKey key(workers, frequency);
    VoltageTableIterator it = table.find(key);
    if(it != table.end()){
        return it->second;
    }else{
        throw runtime_error("Frequency and/or number of virtual cores "
                            "not found in voltage table.");
    }
}

/**
 * Computes a proportional value for power consumption of a configuration.
 * ATTENTION: Assumes linear mapping of the workers on the cores.
 */
static void getPowerProportions(VoltageTable table, uint physicalCores,
                                Frequency frequency, uint coresPerDomain,
                                double& staticPower, double& dynamicPower){
    int numCores = physicalCores;
    double voltage = 0;
    staticPower = 0;
    dynamicPower = 0;
    uint currentCores = 0;
    while(numCores > 0){
        currentCores = (numCores < (int) coresPerDomain)?(uint)numCores:coresPerDomain;
        voltage = getVoltage(table, currentCores, frequency);

        staticPower += voltage;
        dynamicPower += currentCores*frequency*voltage*voltage;

        numCores -= coresPerDomain;
    }
}

RegressionData::RegressionData(const Parameters& p,
                               const FarmConfiguration& configuration,
                               const Smoother<MonitoredSample>* samples):
        _p(p), _configuration(configuration), _samples(samples){
    _topology = _p.mammut.getInstanceTopology();
    _cpus = _topology->getCpus().size();
    _phyCores = _topology->getPhysicalCores().size();
    _phyCoresPerCpu = _topology->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = _topology->getPhysicalCore(0)->getVirtualCores().size();
}

double RegressionData::getUsedPhysicalCores(double numWorkers, bool includeServiceNodes){
    uint numNodes = numWorkers;
    if(includeServiceNodes){
        numNodes += _configuration.getNumServiceNodes();
    }
    return std::min(numNodes, _phyCores);
}

void RegressionData::init(){init(_configuration.getRealValues());}

void RegressionDataServiceTime::init(const KnobsValues& values){
    _numPredictors = 0;
    uint workersCores = values[KNOB_TYPE_WORKERS];
    double physicalCores = getUsedPhysicalCores(workersCores, false);

    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        double frequency = values[KNOB_TYPE_FREQUENCY];
        _invScalFactorFreq = (double)_minFrequency /
                             frequency;
        ++_numPredictors;

        _invScalFactorFreqAndCores = (double)_minFrequency /
                                     (physicalCores * frequency);
        ++_numPredictors;
    }
}

RegressionDataServiceTime::RegressionDataServiceTime(const Parameters& p,
                                                     const FarmConfiguration& configuration,
                                                     const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples),
        _invScalFactorFreq(0),
        _invScalFactorFreqAndCores(0),
        _numPredictors(0){
    _phyCores = _p.mammut.getInstanceTopology()->getPhysicalCores().size();
    _minFrequency = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->
                                                    getAvailableFrequencies().at(0);
    init(_configuration.getRealValues());
}

uint RegressionDataServiceTime::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataServiceTime::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        matrix(rowId++, columnId) = _invScalFactorFreq;
        matrix(rowId++, columnId) = _invScalFactorFreqAndCores;
    }
}

void RegressionDataPower::init(const KnobsValues& values){
    _numPredictors = 0;
    Frequency frequency = values[KNOB_TYPE_FREQUENCY];
    uint workersCores = values[KNOB_TYPE_WORKERS];
    double usedPhysicalCores = getUsedPhysicalCores(workersCores, true);

    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        uint usedCpus = 0;
        if(_p.knobMapping == KNOB_MAPPING_LINEAR){
            switch(_p.knobHyperthreading){
                case KNOB_HT_AUTO:{
                    throw std::runtime_error("This should never happen.");
                }
                case KNOB_HT_LATER:
                case KNOB_HT_SOONER:
                case KNOB_HT_NO:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         (double) _phyCoresPerCpu);
                }break;
                /*
                case KNOB_HT_LATER:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         (double) _phyCoresPerCpu);
                    if(usedPhysicalCores > _phyCores){
                        _additionalContextes = usedPhysicalCores - _phyCores;
                        _additionalContextes = std::min(_additionalContextes, _virtCoresPerPhyCores*_phyCores);
                    }else{
                        _additionalContextes = 0;
                    }
                    ++_numPredictors;
                }break;
                case KNOB_HT_SOONER:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         ((double) _phyCoresPerCpu * (double) _virtCoresPerPhyCores));
                }break;
                */
            }
        }else{
            ; //TODO
        }
        usedCpus = std::min(usedCpus, _cpus);
        uint unusedCpus = _cpus - usedCpus;

        double staticPowerProp = 0, dynamicPowerProp = 0;
        //TODO: Potrebbe non esserci una frequenza se FREQUENCY_NO. In tal caso non possiamo predirre nulla.
        getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                frequency, _phyCoresPerCpu,
                staticPowerProp, dynamicPowerProp);
        _voltagePerUsedSockets = staticPowerProp;
        ++_numPredictors;

        if(_cpus > 1){
            Frequency frequencyUnused;
            switch(_strategyUnused){
                /**
                 * TODO:
                 * For simplicity we assume that the frequency of the unused
                 * sockets is equal to the frequency of the used sockets. This
                 * should be modified as future work.
                 */
                default:{
                    frequencyUnused = frequency;
                }
            }
            double voltage = getVoltage(_p.archData.voltageTable, 0, frequencyUnused);
            _voltagePerUnusedSockets = voltage * unusedCpus;
            ++_numPredictors;
        }

        _dynamicPowerModel = dynamicPowerProp;
        ++_numPredictors;
    }else{
        //TODO: No, moltiplicare comunque per frequenza corrente e v^2, se possibile
        _dynamicPowerModel = usedPhysicalCores;
        ++_numPredictors;
    }
}

RegressionDataPower::RegressionDataPower(const Parameters& p,
                                         const FarmConfiguration& configuration,
                                         const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples), _dynamicPowerModel(0),
        _voltagePerUsedSockets(0), _voltagePerUnusedSockets(0),
        _additionalContextes(0), _numPredictors(0){
    _topology = _p.mammut.getInstanceTopology();
    _cpus = _topology->getCpus().size();
    _phyCores = _topology->getPhysicalCores().size();
    _phyCoresPerCpu = _topology->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = _topology->getPhysicalCore(0)->getVirtualCores().size();
    _strategyUnused = _p.strategyUnusedVirtualCores;
    init(_configuration.getRealValues());
}

uint RegressionDataPower::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataPower::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    matrix(rowId++, columnId) = _dynamicPowerModel;
    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        matrix(rowId++, columnId) = _voltagePerUsedSockets;
        if(_cpus > 1){
            matrix(rowId++, columnId) = _voltagePerUnusedSockets;
        }
    }
    /*
    if(_p.knobHyperthreading != KNOB_HT_NO){
        matrix(rowId++, columnId) = _additionalContextes;
    }
    */
}

Predictor::Predictor(PredictorType type,
                     const Parameters& p,
                     const FarmConfiguration& configuration,
                     const Smoother<MonitoredSample>* samples):
        _type(type), _p(p), _configuration(configuration),
        _samples(samples){
    ;
}

Predictor::~Predictor(){
    ;
}

PredictorLinearRegression::PredictorLinearRegression(PredictorType type,
                                                     const Parameters& p,
                                                     const FarmConfiguration& configuration,
                                                     const Smoother<MonitoredSample>* samples):
        Predictor(type, p, configuration, samples){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            _predictionInput = new RegressionDataServiceTime(p, configuration, samples);
        }break;
        case PREDICTION_POWER:{
            _predictionInput = new RegressionDataPower(p, configuration, samples);
        }break;
    }
    clear();
}

PredictorLinearRegression::~PredictorLinearRegression(){
    clear();
    delete _predictionInput;
}

void PredictorLinearRegression::clear(){
    for(obs_it iterator = _observations.begin();
               iterator != _observations.end();
               iterator++){
        delete iterator->second.data;
    }
    _observations.clear();
    _agingVector.clear();
    _agingVector.reserve(_p.regressionAging);
    _currentAgingId = 0;
}

uint PredictorLinearRegression::getMinimumPointsNeeded(){
    _predictionInput->init(_configuration.getRealValues());
    return std::max<uint>(_predictionInput->getNumPredictors(), 2);
}

double PredictorLinearRegression::getCurrentResponse() const{
    double r = 0.0;
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            r = 1.0 / _samples->average().bandwidth;
        }break;
        case PREDICTION_POWER:{
            r = _samples->average().watts;
        }break;
    }
    return r;
}

bool PredictorLinearRegression::refine(){
    KnobsValues currentValues = _configuration.getRealValues();
    if(_type == PREDICTION_POWER){
        // Add service nodes
        currentValues[KNOB_TYPE_WORKERS] = currentValues[KNOB_TYPE_WORKERS] + _configuration.getNumServiceNodes();
    }
    obs_it lb = _observations.lower_bound(currentValues);

    if(_p.regressionAging){
        if(_agingVector.size() < _p.regressionAging){
            _agingVector.push_back(currentValues);
        }else{
            _agingVector.at(_currentAgingId) = currentValues;
        }
        _currentAgingId = (_currentAgingId + 1) % _p.regressionAging;
    }
    DEBUG("Refining with configuration " << currentValues << ": "
                                         << getCurrentResponse());
    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentValues, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentValues);
        lb->second.data->init();
        lb->second.response = getCurrentResponse();
        return false;
    }else{
        // The key does not exist in the map
        Observation o;
        switch(_type){
            case PREDICTION_BANDWIDTH:{
                o.data = new RegressionDataServiceTime(_p, _configuration, _samples);
            }break;
            case PREDICTION_POWER:{
                o.data = new RegressionDataPower(_p, _configuration, _samples);
            }break;
        }
        o.response = getCurrentResponse();
        _observations.insert(lb, Observations::value_type(currentValues, o));
        return true;
    }
}

void PredictorLinearRegression::prepareForPredictions(){
    if(_observations.size() < getMinimumPointsNeeded()){
        throw std::runtime_error("prepareForPredictions: Not enough " 
                                 "points are present");
    }

    // One observation per column.
    arma::mat dataMl(_observations.begin()->second.data->getNumPredictors(),
                     _observations.size());
    arma::vec responsesMl(_observations.size());

    size_t i = 0;
    for(obs_it iterator = _observations.begin();
               iterator != _observations.end();
               iterator++){
        const Observation& obs = iterator->second;

        if(!_p.regressionAging || contains(_agingVector, iterator->first)){
            obs.data->toArmaRow(i, dataMl);
            responsesMl(i) = obs.response;
            ++i;
        }
    }

    _lr = LinearRegression(dataMl, responsesMl);
    DEBUG("Error in model: " << _lr.ComputeError(dataMl, responsesMl));
    DEBUG("========================== Preparing for predictions: ========================== ");
}

double PredictorLinearRegression::predict(const KnobsValues& values){
    _predictionInput->init(values);

    // One observation per column.
    arma::mat predictionInputMl(_predictionInput->getNumPredictors(), 1);
    arma::vec result(1);
    _predictionInput->toArmaRow(0, predictionInputMl);

    _lr.Predict(predictionInputMl, result);
    if(_type == PREDICTION_BANDWIDTH){
        result.at(0) = 1.0 / result.at(0);
    }

    return result.at(0);
}
    
/**************** PredictorSimple ****************/

PredictorSimple::PredictorSimple(PredictorType type,
                                 const Parameters& p,
                                 const FarmConfiguration& configuration,
                                 const Smoother<MonitoredSample>* samples):
            Predictor(type, p, configuration, samples) {
    Topology* t = _p.mammut.getInstanceTopology();
    _phyCores = t->getPhysicalCores().size();
    _phyCoresPerCpu = t->getCpu(0)->getPhysicalCores().size();
}

double PredictorSimple::getScalingFactor(const KnobsValues& values){
    double usedPhysicalCores = values[KNOB_TYPE_WORKERS];
    return (double)(values[KNOB_TYPE_FREQUENCY] * usedPhysicalCores) /
           (double)(_configuration.getRealValue(KNOB_TYPE_FREQUENCY) *
                    usedPhysicalCores);
}

double PredictorSimple::getPowerPrediction(const KnobsValues& values){
    double staticPower = 0, dynamicPower = 0;
    double usedPhysicalCores = values[KNOB_TYPE_WORKERS];
    getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                        values[KNOB_TYPE_FREQUENCY], _phyCoresPerCpu,
                        staticPower, dynamicPower);
    return dynamicPower;
}

void PredictorSimple::prepareForPredictions(){
    ;
}

double PredictorSimple::predict(const KnobsValues& values){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _samples->average().bandwidth * getScalingFactor(values);
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(values);
        }break;
    }
    return 0.0;
}

Calibrator::Calibrator(const Parameters& p,
                       const FarmConfiguration& configuration,
                       const Smoother<MonitoredSample>* samples):
        _p(p),
        _configuration(configuration),
        _samples(samples),
        _numCalibrationPoints(0),
        _state(CALIBRATION_SEEDS),
        _calibrationStartMs(0),
        _calibrationStartTasks(0),
        _firstPointGenerated(false),
        _primaryPrediction(0), _secondaryPrediction(0),
        _thisPrimary(0), _thisSecondary(0), _noFeasible(false),
        _contractViolations(0), _conservativeValue(_p.conservativeValue){

    PredictorType primary, secondary;
    switch(p.contractType){
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

    switch(p.strategyPrediction){
        case STRATEGY_PREDICTION_SIMPLE:{
            _primaryPredictor = new PredictorSimple(primary, _p, _configuration, _samples);
            _secondaryPredictor = new PredictorSimple(secondary, _p, _configuration, _samples);
        }break;
        case STRATEGY_PREDICTION_REGRESSION_LINEAR:{
            _primaryPredictor = new PredictorLinearRegression(primary, _p, _configuration, _samples);
            _secondaryPredictor = new PredictorLinearRegression(secondary, _p, _configuration, _samples);
        }break;
        default:{
            _primaryPredictor = NULL;
            _secondaryPredictor = NULL;
        }break;
    }

    if(_primaryPredictor && _secondaryPredictor){
        _minNumPoints = std::max(_primaryPredictor->getMinimumPointsNeeded(),
                                 _secondaryPredictor->getMinimumPointsNeeded());
        DEBUG("Minimum number of points required for calibration: " << _minNumPoints);
    }

    _joulesCounter = _localMammut.getInstanceEnergy()->getCounter();
    //TODO Fare meglio con mammut
    //TODO Assicurarsi che il numero totale di configurazioni possibili sia maggiore del numero minimo di punti
}

bool Calibrator::refine(){
    bool p = _primaryPredictor->refine();
    bool s = _secondaryPredictor->refine();
    return p & s;
}

bool Calibrator::isAccurate(double primaryValue,
                                      double secondaryValue) const{
    double primaryError = std::abs((primaryValue - _primaryPrediction)/
                                   primaryValue)*100.0;
    double secondaryError = 0;
    secondaryError = std::abs((secondaryValue - _secondaryPrediction)/
                               secondaryValue)*100.0;

    DEBUG("Primary prediction: " << _primaryPrediction << " " <<
          "Secondary prediction: " << _secondaryPrediction);
    DEBUG("Primary error: " << primaryError << " " <<
          "Secondary error: " << secondaryError);

    if(primaryError > _p.maxPrimaryPredictionError ||
       secondaryError > _p.maxSecondaryPredictionError){
        return false;
    }else{
        return true;
    }
}

bool Calibrator::isBestSuboptimalValue(double x, double y) const{
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

bool Calibrator::isBestSecondaryValue(double x, double y) const{
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

bool Calibrator::isFeasiblePrimaryValue(double value, bool conservative) const{
    double conservativeOffset = 0;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            if(conservative){
                conservativeOffset = ((_p.overloadThresholdFarm - _p.underloadThresholdFarm) *
                                       _conservativeValue) / 100.0;
            }
            return value > _p.underloadThresholdFarm + conservativeOffset &&
                   value < _p.overloadThresholdFarm - conservativeOffset;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            if(conservative){
                conservativeOffset = (_p.requiredBandwidth * _conservativeValue) / 100.0;
            }
            return value > _p.requiredBandwidth + conservativeOffset;
        }break;
        case CONTRACT_POWER_BUDGET:{
            if(conservative){
                conservativeOffset = (_p.powerBudget * _conservativeValue) / 100.0;
            }
            return value < _p.powerBudget - conservativeOffset;
        }break;
        default:{
            return false;
        }break;
    }
    return false;
}

KnobsValues Calibrator::getBestKnobsValues(double primaryValue){
    KnobsValues bestValues(KNOB_VALUE_REAL);
    KnobsValues bestSuboptimalValues = _configuration.getRealValues();

    double primaryPrediction = 0;
    double secondaryPrediction = 0;

    double bestPrimaryPrediction = 0;
    double bestSecondaryPrediction = 0;
    double bestSuboptimalValue = primaryValue;

    bool feasibleSolutionFound = false;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            // We have to minimize the power.
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

    vector<KnobsValues> combinations = _configuration.getAllRealCombinations();
    for(size_t i = 0; i < combinations.size(); i++){
        KnobsValues currentValues = combinations.at(i);
        primaryPrediction = _primaryPredictor->predict(currentValues);
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                primaryPrediction = (_samples->average().bandwidth / primaryPrediction) *
                                     _samples->average().utilization;
            }break;
            default:{
                ;
            }
        }

        //std::cout << currentValues << " " << primaryPrediction << " ";
        if(isFeasiblePrimaryValue(primaryPrediction, true)){
            secondaryPrediction = _secondaryPredictor->predict(currentValues);
            //std::cout << secondaryPrediction;
            if(isBestSecondaryValue(secondaryPrediction, bestSecondaryPrediction)){
                bestValues = currentValues;
                feasibleSolutionFound = true;
                bestPrimaryPrediction = primaryPrediction;
                bestSecondaryPrediction = secondaryPrediction;
            }
        }else if(!feasibleSolutionFound &&
                 isBestSuboptimalValue(primaryPrediction, bestSuboptimalValue)){
            bestSuboptimalValue = primaryPrediction;
            bestSuboptimalValues = currentValues;
        }
        //std::cout << std::endl;
    }

    if(feasibleSolutionFound){
        DEBUG("Best solution found " << bestValues);
        _primaryPrediction = bestPrimaryPrediction;
        _secondaryPrediction = bestSecondaryPrediction;
        DEBUG("Primary prediction: " << _primaryPrediction);
        DEBUG("Secondary prediction: " << _secondaryPrediction);
        return bestValues;
    }else{
        DEBUG("Suboptimal solution found.");
        _primaryPrediction = bestSuboptimalValue;
        _secondaryPrediction = NO_PREDICTION;
        return bestSuboptimalValues;
    }
}

bool Calibrator::isContractViolated(double primaryValue) const{
    return !isFeasiblePrimaryValue(primaryValue, false);
}

bool Calibrator::phaseChanged(double primaryValue, double secondaryValue) const{
    return (_thisPrimary != -1 && (primaryValue > 2*_thisPrimary || primaryValue < 0.5*_thisPrimary)) ||
        (_thisSecondary != -1 && (secondaryValue > 2*_thisSecondary || secondaryValue < 0.5*_thisSecondary));
}

void Calibrator::startCalibrationStat(uint64_t totalTasks){
    _numCalibrationPoints = 0;
    _calibrationStartMs = getMillisecondsTime();
    _calibrationStartTasks = totalTasks;
    if(_joulesCounter){
        _joulesCounter->reset();
    }
}

void Calibrator::stopCalibrationStat(uint64_t totalTasks){
    if(_numCalibrationPoints){
        CalibrationStats cs;
        cs.numSteps = _numCalibrationPoints;
        cs.duration = (getMillisecondsTime() - _calibrationStartMs);
        cs.numTasks = totalTasks - _calibrationStartTasks;
        if(_joulesCounter){
            cs.joules = ((CounterCpus*) _joulesCounter)->getJoulesCoresAll();
        }
        _calibrationStats.push_back(cs);
        _numCalibrationPoints = 0;
    }
}

void Calibrator::updatePredictions(const KnobsValues& next){
    KnobsValues real;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        double realv;
        _configuration.getKnob((KnobType)i)->getRealFromRelative(next[(KnobType)i], realv);
        real[(KnobType)i] = realv;
    }
    _primaryPredictor->prepareForPredictions();
    _secondaryPredictor->prepareForPredictions();
    _primaryPrediction = _primaryPredictor->predict(real);
    _secondaryPrediction = _secondaryPredictor->predict(real);
}

void Calibrator::updateConservativeValue(){
    /*
    if(_contractViolations >= 1){
        _conservativeValue = _p.conservativeValue + _contractViolations;
    }*/
    switch(_p.contractType){
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_UTILIZATION:{
            _conservativeValue = _samples->coefficientVariation().bandwidth;
        }break;
        case CONTRACT_POWER_BUDGET:{
            _conservativeValue = _samples->coefficientVariation().watts;
        }break;
    }
}

KnobsValues Calibrator::getNextKnobsValues(double primaryValue,
                                           double secondaryValue,
                                           u_int64_t totalTasks){
    KnobsValues kv;
    bool contractViolated = isContractViolated(primaryValue);

    /**
     * The first point is generated as soon as the application starts.
     * Accordingly, we do not executed tasks in the original configuration
     * used to create the application. For this reason, we do not use
     * it to refine the model.
     * E.g. The application has been created with configuration X
     * and as soon as it starts we move it to configuration Y.
     * We do not refine with configuration X since it has never
     * been real executed.
     **/
    if(_firstPointGenerated){
        /**
         * When we are in the final state and the contract is not violated
         * (so we remain in the same configuration), we do not still have to
         * refine the predictions.
         */
        if(_state == CALIBRATION_FINISHED){
            ;
        }else{
            refine();
            ++_numCalibrationPoints;
        }
    }else{
        _firstPointGenerated = true;
        startCalibrationStat(totalTasks);
    }

    switch(_state){
        case CALIBRATION_SEEDS:{
            _noFeasible = false;
            if(_numCalibrationPoints < _minNumPoints){
                kv = generateRelativeKnobsValues();
            }else{
                if(!isAccurate(primaryValue, secondaryValue) ||
                   _numCalibrationPoints == _minNumPoints){
                    kv = generateRelativeKnobsValues();
                    updatePredictions(kv);
                    DEBUG("[Calibrator]: High prediction error. Adding new seed.");
                }else{
                    kv = getBestKnobsValues(primaryValue);
                    DEBUG("[Calibrator]: Low error but the contract is violated, adjusting.");

                    if(_secondaryPrediction == NO_PREDICTION){
                        DEBUG("No feasible solutions found.");
                        _noFeasible = true;
                    }

                    _state = CALIBRATION_FINISHED;
                    _thisPrimary = -1; // primaryValue;
                    _thisSecondary = -1; //secondaryValue;

                    DEBUG("[Calibrator]: Moving to finished");
                    DEBUG("[Calibrator]: Finished in " << _numCalibrationPoints <<
                          " steps with configuration " << kv);
                    stopCalibrationStat(totalTasks);
                }
            }
        }break;
        case CALIBRATION_FINISHED:{
            updateConservativeValue();

            if(phaseChanged(primaryValue, secondaryValue)){
                kv = reset();

                _state = CALIBRATION_SEEDS;
                startCalibrationStat(totalTasks);

                _contractViolations = 0;
            }else if((!_noFeasible && contractViolated) ||
                     !isAccurate(primaryValue, secondaryValue)){
                kv = generateRelativeKnobsValues();
                updatePredictions(kv);

                refine();

                _state = CALIBRATION_SEEDS;
                startCalibrationStat(totalTasks);

                ++_contractViolations;
            }else{
                kv = _configuration.getRealValues();
            }
        }break;
    }

    return kv;
}

std::vector<CalibrationStats> Calibrator::getCalibrationsStats() const{
    return _calibrationStats;
}

bool Calibrator::isCalibrating() const{
    return _state != CALIBRATION_FINISHED;
}

KnobsValues CalibratorDummy::getNextKnobsValues(double primaryValue,
                                                double secondaryValue,
                                                u_int64_t totalTasks){
    return getBestKnobsValues(primaryValue);
}


CalibratorRandom::CalibratorRandom(const Parameters& p,
                    const FarmConfiguration& configuration,
                    const Smoother<MonitoredSample>* samples):
                        Calibrator(p, configuration, samples){
    srand(time(NULL));
}

KnobsValues CalibratorRandom::generateRelativeKnobsValues() const{
    KnobsValues r(KNOB_VALUE_RELATIVE);
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(_configuration.getKnob((KnobType) i)->needsCalibration()){
            r[(KnobType)i] = rand() % 100;
        }else{
            /**
             * If we do not need to automatically find the value for this knob,
             * then it has only 0 or 1 possible value. Accordingly, we can set
             * it to any value. Here we set it to 100.0 for readability.
             */
            r[(KnobType)i] = 100.0;
        }
    }
    return r;
}


CalibratorLowDiscrepancy::CalibratorLowDiscrepancy(const Parameters& p,
                                                   const FarmConfiguration& configuration,
                                                   const Smoother<MonitoredSample>* samples):
        Calibrator(p, configuration, samples){
    uint d = 0;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(_configuration.getKnob((KnobType) i)->needsCalibration()){
            ++d;
        }
    }
    DEBUG("[Calibrator] Will generate low discrepancy points in " << d <<
          " dimensions");
    const gsl_qrng_type* generatorType;
    switch(_p.strategyCalibration){
        case STRATEGY_CALIBRATION_NIEDERREITER:{
            generatorType = gsl_qrng_niederreiter_2;
        }break;
        case STRATEGY_CALIBRATION_SOBOL:{
            generatorType = gsl_qrng_sobol;
        }break;
        case STRATEGY_CALIBRATION_HALTON:{
            generatorType = gsl_qrng_halton;
        }break;
        case STRATEGY_CALIBRATION_HALTON_REVERSE:{
            generatorType = gsl_qrng_reversehalton;
        }break;
        default:{
            throw std::runtime_error("CalibratorLowDiscrepancy: Unknown "
                                     "generator type: " +
                                     _p.strategyCalibration);
        }break;
    }
    _generator = gsl_qrng_alloc(generatorType, d);
    _normalizedPoint = new double[d];
}

CalibratorLowDiscrepancy::~CalibratorLowDiscrepancy(){
    gsl_qrng_free(_generator);
    delete[] _normalizedPoint;
}

KnobsValues CalibratorLowDiscrepancy::generateRelativeKnobsValues() const{
    KnobsValues r(KNOB_VALUE_RELATIVE);
    gsl_qrng_get(_generator, _normalizedPoint);
    size_t nextCoordinate = 0;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(_configuration.getKnob((KnobType) i)->needsCalibration()){
            r[(KnobType)i] = _normalizedPoint[nextCoordinate]*100.0;
            ++nextCoordinate;
        }else{
            /**
             * If we do not need to automatically find the value for this knob,
             * then it has only 0 or 1 possible value. Accordingly, we can set
             * it to any value. Here we set it to 100.0 for readability.
             */
            r[(KnobType)i] = 100.0;
        }
    }
    return r;
}

KnobsValues CalibratorLowDiscrepancy::reset(){
    gsl_qrng_init(_generator);
    KnobsValues kv = generateRelativeKnobsValues();
    _primaryPredictor->clear();
    _secondaryPredictor->clear();
    DEBUG("[Calibrator]: Moving to seeds");
    return kv;
}

Frequency CalibratorLiMartinez::findNearestFrequency(Frequency f) const{
    Frequency bestDistance = _availableFrequencies.back();
    Frequency bestFrequency = _availableFrequencies.back();
    for(size_t i = 0; i < _availableFrequencies.size(); i++){
        Frequency distance = std::abs(_availableFrequencies.at(i) - f);
        if(distance < bestDistance){
            bestDistance = distance;
            bestFrequency = _availableFrequencies.at(i);
        }
    }
    return bestFrequency;
}

void CalibratorLiMartinez::goRight(){
    _low1 = _low2;
    _high1 = _mid2 - 1;
    _low2 = _mid2 + 1;
    _mid1 = (_low1 + _high1) / 2.0;
    _mid2 = (_low2 + _high2) / 2.0;
}

void CalibratorLiMartinez::goLeft(){
    _high1 = _mid1 - 1;
    _low2 = _mid1 + 1;
    _high2 = _high1;
    _mid1 = (_low1 + _high1) / 2.0;
    _mid2 = (_low2 + _high2) / 2.0;
}

CalibratorLiMartinez::CalibratorLiMartinez(const Parameters& p,
                                           const FarmConfiguration& configuration,
                                           const Smoother<MonitoredSample>* samples):
    Calibrator(p, configuration, samples),
    _firstPointGenerated(false), _low1(0), _mid1(0), _high1(0),
    _low2(0), _mid2(0), _high2(0), _midId(1),
    _availableFrequencies(_p.mammut.getInstanceCpuFreq()->getDomains().back()->getAvailableFrequencies()),
    _currentWatts(DBL_MAX), _optimalWatts(DBL_MAX),
    _optimalFrequency(_availableFrequencies.back()), _optimalWorkers(1),
    _currentBw(0),_leftBw(0), _rightBw(0), _optimalFound(false), _improved(false){
    ;
}

CalibratorLiMartinez::~CalibratorLiMartinez(){
    ;
}

KnobsValues CalibratorLiMartinez::getNextKnobsValues(double primaryValue,
                                                     double secondaryValue,
                                                     u_int64_t totalTasks){
    KnobsValues kv(KNOB_VALUE_REAL);
    kv[KNOB_TYPE_MAPPING] = KNOB_MAPPING_LINEAR;

    if(!_firstPointGenerated){
        _firstPointGenerated = true;
        uint maxWorkers = _configuration.getKnob(KNOB_TYPE_WORKERS)->getRealValue();
        _low2 = 1;
        _mid2 = maxWorkers / 2.0;
        _high2 = maxWorkers;

        kv[KNOB_TYPE_WORKERS] = _mid2;
        kv[KNOB_TYPE_FREQUENCY] = _availableFrequencies.back();
        _midId = 2;

        startCalibrationStat(totalTasks);

        DEBUG("Generating first point: " << kv);
        ++_numCalibrationPoints;
    }else{
        if(_optimalFound){
            return _optimalKv;
        }else if(!isContractViolated(primaryValue)){
            ++_numCalibrationPoints;
            _currentWatts = secondaryValue;

            if(_currentWatts < _optimalWatts){
                _improved = true;
                DEBUG("Found a new optimal watts: " << _currentWatts << " vs. " << _optimalWatts);
                _optimalWatts = _currentWatts;
                _optimalFrequency = _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealValue();
                _optimalWorkers = _configuration.getKnob(KNOB_TYPE_WORKERS)->getRealValue();
                DEBUG("Optimal: " << _optimalWorkers << ", " << _optimalFrequency);
            }

            // We should keep decreasing the frequency
            Frequency currentFrequency = _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealValue();
            kv[KNOB_TYPE_WORKERS] = _configuration.getKnob(KNOB_TYPE_WORKERS)->getRealValue();

            Frequency nextFrequency = currentFrequency * (_p.requiredBandwidth / primaryValue);
            nextFrequency = findNearestFrequency(nextFrequency);
            DEBUG("Required BW: " << _p.requiredBandwidth << " Current BW: " << primaryValue << " Best frequency: " << nextFrequency);
            if(nextFrequency == currentFrequency){
                --_numCalibrationPoints;
                goto changeworkers;
            }else{
                kv[KNOB_TYPE_FREQUENCY] = nextFrequency;
                DEBUG("Keeping going down on frequencies. We move to: " << kv);
            }
        }else{
changeworkers:
            ++_numCalibrationPoints;
            // I have to change the number of workers
            kv[KNOB_TYPE_FREQUENCY] = _availableFrequencies.back();

            if(_optimalWatts == DBL_MAX){
                // Still I have not found a number of workers that satisfied
                // the time requirement. I increase workers. (Go right).
                kv[KNOB_TYPE_WORKERS] = _mid2;
                goRight();
                _midId = 2;
            }else if(_currentWatts > _optimalWatts || !_improved){
                DEBUG("This number of workers is worst than the best we found "
                      "up to now.");
                // This number of workers is not ok
                if(_midId == 1){
                    kv[KNOB_TYPE_WORKERS] = _mid2;
                    _midId = 2;
                    DEBUG("Trying with the right side. We move to " << kv);
                }else{
                    // Both explored and both are not ok, finished
                    kv[KNOB_TYPE_WORKERS] = _optimalWorkers;
                    kv[KNOB_TYPE_FREQUENCY] = _optimalFrequency;
                    _optimalFound = true;
                    _optimalKv = kv;
                    stopCalibrationStat(totalTasks);
                    DEBUG("Both side are worst. Terminated with: " << kv);
                }
            }else{
                _improved = false;
                if(_midId == 1){
                    goLeft();
                }else{
                    goRight();
                }

                DEBUG("New interval 1: [" << _low1 << "," << _mid1 << "," << _high1 << "]");
                DEBUG("New interval 2: [" << _low2 << "," << _mid2 << "," << _high2 << "]");

                if(_low1 <= _high1){
                    _midId = 1;
                    kv[KNOB_TYPE_WORKERS] = _mid1;
                    DEBUG("We move to " << kv);
                }else if(_low2 <= _high2){
                    _midId = 2;
                    kv[KNOB_TYPE_WORKERS] = _mid2;
                    DEBUG("We move to " << kv);
                }else{
                    kv[KNOB_TYPE_WORKERS] = _optimalWorkers;
                    kv[KNOB_TYPE_FREQUENCY] = _optimalFrequency;
                    _optimalFound = true;
                    _optimalKv = kv;
                    stopCalibrationStat(totalTasks);
                    DEBUG("Exploration finished with: " << kv);
                }
            }
        }
    }
    return kv;
}

}

