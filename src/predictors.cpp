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


namespace adpff{

using namespace mammut::cpufreq;
using namespace mammut::utils;

static double getVoltage(VoltageTable table, const KnobsValues& values){
    VoltageTableKey key(values[KNOB_TYPE_WORKERS], values[KNOB_TYPE_FREQUENCY]);
    VoltageTableIterator it = table.find(key);
    if(it != table.end()){
        return it->second;
    }else{
        throw runtime_error("Frequency and/or number of virtual cores "
                                 "not found in voltage table.");
    }
}

RegressionData::RegressionData(const Parameters& p,
                               const FarmConfiguration& configuration,
                               const Smoother<MonitoredSample>* samples):
        _p(p), _configuration(configuration), _samples(samples){
    ;
}

void RegressionData::init(){init(_configuration.getRealValues());}

void RegressionDataServiceTime::init(const KnobsValues& values){
    _numPredictors = 0;
    double physicalCores = std::min((uint) values[KNOB_TYPE_WORKERS], _phyCores);

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
    uint numWorkers = values[KNOB_TYPE_WORKERS];
    Frequency frequency = values[KNOB_TYPE_FREQUENCY];
    double usedPhysicalCores = std::min(numWorkers, _phyCores);

    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        double voltage = getVoltage(_p.archData.voltageTable, values);
        uint usedCpus = 0;
        if(_p.knobMapping == KNOB_MAPPING_LINEAR){
            switch(_p.knobHyperthreading){
                case KNOB_HT_AUTO:{
                    throw std::runtime_error("This should never happen.");
                }
                case KNOB_HT_NO:{
                    usedCpus = std::ceil((double) numWorkers /
                                         (double) _phyCoresPerCpu);
                }break;
                case KNOB_HT_LATER:{
                    usedCpus = std::ceil((double) numWorkers /
                                         (double) _phyCoresPerCpu);
                    if(numWorkers > _phyCores){
                        _additionalContextes = numWorkers - _phyCores;
                        _additionalContextes = std::min(_additionalContextes, _virtCoresPerPhyCores*_phyCores);
                    }else{
                        _additionalContextes = 0;
                    }
                    ++_numPredictors;
                }break;
                case KNOB_HT_SOONER:{
                    usedCpus = std::ceil((double) numWorkers /
                                         ((double) _phyCoresPerCpu * (double) _virtCoresPerPhyCores));
                }break;
            }
        }else{
            ; //TODO
        }
        usedCpus = std::min(usedCpus, _cpus);
        uint unusedCpus = _cpus - usedCpus;

        _voltagePerUsedSockets = voltage * (double) usedCpus;
        ++_numPredictors;

        if(_cpus > 1){
            _voltagePerUnusedSockets = voltage * (double) unusedCpus;
            ++_numPredictors;
        }

        _dynamicPowerModel = (usedPhysicalCores*frequency*voltage*voltage);
        ++_numPredictors;
    }else{
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
    Topology* t = _p.mammut.getInstanceTopology();
    _cpus = t->getCpus().size();
    _phyCores = t->getPhysicalCores().size();
    _phyCoresPerCpu = t->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = t->getPhysicalCore(0)->getVirtualCores().size();
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
    if(_p.knobHyperthreading != KNOB_HT_NO){
        matrix(rowId++, columnId) = _additionalContextes;
    }
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

void PredictorLinearRegression::refine(){
    const KnobsValues& currentValues = _configuration.getRealValues();
    obs_it lb = _observations.lower_bound(currentValues);

    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentValues, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentValues);
        lb->second.data->init();
        lb->second.response = getCurrentResponse();
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
    }

    DEBUG("Refining with configuration " << currentValues << ": "
                                         << getCurrentResponse());
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
        obs.data->toArmaRow(i, dataMl);
        responsesMl(i) = obs.response;
        ++i;
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
    ;
}

double PredictorSimple::getScalingFactor(const KnobsValues& values){
    return (double)(values[KNOB_TYPE_FREQUENCY] * values[KNOB_TYPE_WORKERS]) /
           (double)(_configuration.getRealValue(KNOB_TYPE_FREQUENCY) *
                    _configuration.getRealValue(KNOB_TYPE_WORKERS));
}

double PredictorSimple::getPowerPrediction(const KnobsValues& values){
    Voltage v = getVoltage(_p.archData.voltageTable, values);
    return values[KNOB_TYPE_WORKERS]*values[KNOB_TYPE_FREQUENCY]*v*v;
}

void PredictorSimple::prepareForPredictions(){
    ;
}

double PredictorSimple::predict(const KnobsValues& values){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _samples->average().bandwidth *
                   getScalingFactor(values);
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
        _state(CALIBRATION_SEEDS),
        _numCalibrationPoints(0), _calibrationStartMs(0),
        _firstPointGenerated(false), _primaryPrediction(0),
        _secondaryPrediction(0){

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
            ;
        }break;
    }

    _minNumPoints = std::max(_primaryPredictor->getMinimumPointsNeeded(),
                             _secondaryPredictor->getMinimumPointsNeeded());
    DEBUG("Minimum number of points required for calibration: " << _minNumPoints);
    //TODO Assicurarsi che il numero totale di configurazioni possibili sia maggiore del numero minimo di punti
}

void Calibrator::refine(bool isContractViolated){
    // We have to refine only if a new configuration will be generated.
    if(_state == CALIBRATION_FINISHED && isContractViolated){
        return;
    }
    _primaryPredictor->refine();
    _secondaryPredictor->refine();
}

bool Calibrator::highError(double primaryValue, double secondaryValue) const{
    double primaryError = std::abs((primaryValue - _primaryPrediction)/
                                   primaryValue)*100.0;
    double secondaryError;

    // TODO: This check now works because both service time and power are always  > 0
    // In the future we must find another way to indicat that secondary prediction
    // has not been done.
    if(_secondaryPrediction >= 0.0){
        secondaryError = std::abs((secondaryValue - _secondaryPrediction)/
                                   secondaryValue)*100.0;
    }else{
        secondaryError = 0.0;
    }
    DEBUG("Primary prediction: " << _primaryPrediction << " " <<
          "Secondary prediction: " << _secondaryPrediction);
    DEBUG("Primary error: " << primaryError << " " <<
          "Secondary error: " << secondaryError);
    return primaryError > _p.maxPrimaryPredictionError ||
           secondaryError > _p.maxSecondaryPredictionError;
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

bool Calibrator::isFeasiblePrimaryValue(double value, double tolerance) const{
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

KnobsValues Calibrator::getBestKnobsValues(double primaryValue,
                                           double secondaryValue,
                                           u_int64_t remainingTasks){
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
                remainingTime = (double) remainingTasks / primaryPrediction;
            }break;
            case CONTRACT_PERF_UTILIZATION:{
                primaryPrediction = (_samples->average().bandwidth / primaryPrediction) *
                                     _samples->average().utilization;
            }break;
            default:{
                ;
            }
        }

        if(isFeasiblePrimaryValue(primaryPrediction)){
            secondaryPrediction = _secondaryPredictor->predict(currentValues);
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
        // In the future we must find another way to indicate that secondary prediction
        // has not been done.
        _secondaryPrediction = -1;
        return bestSuboptimalValues;
    }
}

bool Calibrator::isContractViolated(double primaryValue) const{
    double tolerance = 0;
    double maxError = _p.maxPrimaryPredictionError;
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
    return !isFeasiblePrimaryValue(primaryValue, tolerance);
}

KnobsValues Calibrator::getNextKnobsValues(double primaryValue,
                                           double secondaryValue,
                                           u_int64_t remainingTasks){
    KnobsValues kv;

    if(_numCalibrationPoints == 0){
        _calibrationStartMs = getMillisecondsTime();
    }

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
        refine(isContractViolated(primaryValue));
    }else{
        _firstPointGenerated = true;
    }

    switch(_state){
        case CALIBRATION_SEEDS:{
            kv = generateRelativeKnobsValues();
            if(_numCalibrationPoints + 1 >= _minNumPoints){
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: Moving to predict");
            }
        }break;
        case CALIBRATION_TRY_PREDICT:{
            kv = getBestKnobsValues(primaryValue, secondaryValue,
                                    remainingTasks);
            _state = CALIBRATION_EXTRA_POINT;
            DEBUG("[Calibrator]: Moving to extra");
        }break;
        case CALIBRATION_EXTRA_POINT:{
            if(highError(primaryValue, secondaryValue)){
                kv = generateRelativeKnobsValues();
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: High error");
                DEBUG("[Calibrator]: Moving to predict");
            }else if(isContractViolated(primaryValue)){
                DEBUG("[Calibrator]: Contract violated");
                kv = getBestKnobsValues(primaryValue, secondaryValue,
                                        remainingTasks);
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: Moving to predict");
            }else{
                kv = getBestKnobsValues(primaryValue, secondaryValue,
                                        remainingTasks);
                _state = CALIBRATION_FINISHED;
                _primaryPredictor->clear();
                _secondaryPredictor->clear();

                CalibrationStats cs;
                // We do -1 because we counted the current point and now we
                // discovered it isn't a calibration point.
                cs.numSteps = _numCalibrationPoints - 1;
                cs.duration = (getMillisecondsTime() - _calibrationStartMs);
                _calibrationStats.push_back(cs);
                DEBUG("[Calibrator]: Moving to finished");
                DEBUG("[Calibrator]: Finished in " << _numCalibrationPoints - 1 <<
                      " steps with configuration " << kv);
            }
        }break;
        case CALIBRATION_FINISHED:{
            /*
            if(highError()){
                fc = _seeds.at(0);
                _nextSeed = 1;
                _state = CALIBRATION_SEEDS;
                DEBUG("========Moving to seeds");
            }else*/
            if(isContractViolated(primaryValue)){
                reset();
                kv = generateRelativeKnobsValues();
                _numCalibrationPoints = 1;
                _state = CALIBRATION_SEEDS;
                _calibrationStartMs = getMillisecondsTime();
                DEBUG("[Calibrator]: Moving to seeds");
            }else{
                kv = _configuration.getRealValues();
            }
        }break;
    }

    if(_state != CALIBRATION_FINISHED){
        ++_numCalibrationPoints;
    }

    return kv;
}

std::vector<CalibrationStats> Calibrator::getCalibrationsStats() const{
    return _calibrationStats;
}

CalibratorLowDiscrepancy::CalibratorLowDiscrepancy(const Parameters& p,
                                                   const FarmConfiguration& configuration,
                                                   const Smoother<MonitoredSample>* samples):
        Calibrator(p, configuration, samples){
    uint d = 0;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(_configuration.getKnob((KnobType) i)->autoFind()){
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
        if(_configuration.getKnob((KnobType) i)->autoFind()){
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

void CalibratorLowDiscrepancy::reset(){
    gsl_qrng_init(_generator);
}

}

