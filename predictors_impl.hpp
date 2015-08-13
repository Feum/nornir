/*
 * predictors_impl.hpp
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

/*! \file predictors_impl.hpp
 * \brief Predictors used by adaptive farm.
 *
 */
#ifndef PREDICTORS_IMPL_HPP_
#define PREDICTORS_IMPL_HPP_

#include "predictors.hpp"
#include "manager.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_PREDICTORS
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#define DEBUGB(x) do {x} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif


namespace adpff{

using namespace mammut::cpufreq;
using namespace mammut::utils;

void RegressionData::init(){init(_manager._currentConfiguration);}

void RegressionDataServiceTime::init(const FarmConfiguration& configuration){
    _numPredictors = 0;
    double usedPhysicalCores = std::min(configuration.numWorkers, _manager._numPhysicalCores);

    _physicalCores = usedPhysicalCores;
    ++_numPredictors;

    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        _invScalFactorPhysical = 1.0 / (usedPhysicalCores * ((double)configuration.frequency /
                                                             (double)_manager._availableFrequencies.at(0)));
        ++_numPredictors;
    }

    if(_manager._p.strategyHyperthreading != STRATEGY_HT_NO){
        _workers = (double)configuration.numWorkers;
        ++_numPredictors;

        _invScalFactorWorkers = 1.0 / ((double) configuration.numWorkers * ((double)configuration.frequency /
                                                                            (double)_manager._availableFrequencies.at(0)));
        ++_numPredictors;
    }
    /*
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        _frequency = configuration.frequency;
        ++_numPredictors;
    }
    */
}

RegressionDataServiceTime::RegressionDataServiceTime(const ManagerFarm& manager):
        RegressionData(manager),
        _physicalCores(0), _invScalFactorPhysical(0),
        _workers(0), _invScalFactorWorkers(0),
        _frequency(0), _numPredictors(0){
    init(manager._currentConfiguration);
}

RegressionDataServiceTime::RegressionDataServiceTime(const ManagerFarm& manager,
                                                     const FarmConfiguration& configuration):
        RegressionData(manager){
    init(configuration);
}

uint RegressionDataServiceTime::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataServiceTime::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    //TODO: !!!
    matrix(rowId++, columnId) = _physicalCores;
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        matrix(rowId++, columnId) = _invScalFactorPhysical;
    }
    if(_manager._p.strategyHyperthreading != STRATEGY_HT_NO){
        matrix(rowId++, columnId) = _workers;
        matrix(rowId++, columnId) = _invScalFactorWorkers;
    }
    /*
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        matrix(rowId++, columnId) = _frequency;
    }
    */
}

void RegressionDataPower::init(const FarmConfiguration& configuration){
    _numPredictors = 0;
    double usedPhysicalCores = std::min(configuration.numWorkers, _manager._numPhysicalCores);

    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        double voltage = _manager.getVoltage(configuration);
        uint usedCpus;
        if(_manager._p.strategyMapping == STRATEGY_MAPPING_LINEAR){
            switch(_manager._p.strategyHyperthreading){
                case STRATEGY_HT_NO:{
                    usedCpus = std::ceil((double) configuration.numWorkers /
                                         (double) _manager._numPhysicalCoresPerCpu);
                }break;
                case STRATEGY_HT_YES_LATER:{
                    usedCpus = std::ceil((double) configuration.numWorkers /
                                         (double) _manager._numPhysicalCoresPerCpu);
                    if(configuration.numWorkers > _manager._numPhysicalCores){
                        _additionalContextes = configuration.numWorkers - _manager._numPhysicalCores;
                        _additionalContextes = std::min(_additionalContextes, _manager._numVirtualCoresPerPhysicalCore*_manager._numPhysicalCores);
                    }else{
                        _additionalContextes = 0;
                    }
                    ++_numPredictors;
                }break;
                case STRATEGY_HT_YES_SOONER:{
                    usedCpus = std::ceil((double) configuration.numWorkers /
                                         ((double) _manager._numPhysicalCoresPerCpu * (double) _manager._numVirtualCoresPerPhysicalCore));
                }break;
            }
        }else{
            ; //TODO
        }
        usedCpus = std::min(usedCpus, _manager._numCpus);
        uint unusedCpus = _manager._numCpus - usedCpus;

        _voltagePerUsedSockets = voltage * (double) usedCpus;
        ++_numPredictors;

        if(_manager._numCpus > 1){
            _voltagePerUnusedSockets = voltage * (double) unusedCpus;
            ++_numPredictors;
        }

        _dynamicPowerModel = (usedPhysicalCores*configuration.frequency*voltage*voltage);
        ++_numPredictors;
    }else{
        _dynamicPowerModel = usedPhysicalCores;
        ++_numPredictors;
    }
}

RegressionDataPower::RegressionDataPower(const ManagerFarm& manager):
        RegressionData(manager), _dynamicPowerModel(0),
        _voltagePerUsedSockets(0), _voltagePerUnusedSockets(0),
        _additionalContextes(0), _numPredictors(0){
    init(_manager._currentConfiguration);
}

RegressionDataPower::RegressionDataPower(const ManagerFarm& manager,
                                         const FarmConfiguration& configuration):
        RegressionData(manager){
    init(configuration);
}

uint RegressionDataPower::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataPower::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    matrix(rowId++, columnId) = _dynamicPowerModel;
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        matrix(rowId++, columnId) = _voltagePerUsedSockets;
        if(_manager._numCpus > 1){
            matrix(rowId++, columnId) = _voltagePerUnusedSockets;
        }
    }
    if(_manager._p.strategyHyperthreading != STRATEGY_HT_NO){
        matrix(rowId++, columnId) = _additionalContextes;
    }
}

PredictorLinearRegression::PredictorLinearRegression(PredictorType type,
                                                     const ManagerFarm& manager):
                                                    _type(type), _manager(manager){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            _predictionInput = new RegressionDataServiceTime(_manager);
        }break;
        case PREDICTION_POWER:{
            _predictionInput = new RegressionDataPower(_manager);
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
    //TODO: Brutto a vedersi, ripulire
    _predictionInput->init(_manager._currentConfiguration);
    return std::max<uint>(_predictionInput->getNumPredictors(), 2);
}

double PredictorLinearRegression::getCurrentResponse() const{
    double r = 0.0;
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            r = 1.0 / _manager._samples->average().bandwidth;
        }break;
        case PREDICTION_POWER:{
            r = _manager._samples->average().watts.cores;
        }break;
    }
    return r;
}

void PredictorLinearRegression::refine(){
    const FarmConfiguration& currentConf = _manager._currentConfiguration;
    obs_it lb = _observations.lower_bound(currentConf);

    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentConf, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentConf);
        lb->second.data->init();
        lb->second.response = getCurrentResponse();
    }else{
        // The key does not exist in the map
        Observation o;
        switch(_type){
            case PREDICTION_BANDWIDTH:{
                o.data = new RegressionDataServiceTime(_manager);
            }break;
            case PREDICTION_POWER:{
                o.data = new RegressionDataPower(_manager);
            }break;
        }
        o.response = getCurrentResponse();
        _observations.insert(lb, Observations::value_type(currentConf, o));
    }

    DEBUG("Refining with configuration " << currentConf << ": "
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

double PredictorLinearRegression::predict(const FarmConfiguration& configuration){
    _predictionInput->init(configuration);

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
                                 const ManagerFarm& manager):
    _type(type), _manager(manager){
    ;
}

double PredictorSimple::getScalingFactor(const FarmConfiguration& configuration){
    return (double)(configuration.frequency * configuration.numWorkers) /
           (double)(_manager._currentConfiguration.frequency * _manager._currentConfiguration.numWorkers);
}

double PredictorSimple::getPowerPrediction(const FarmConfiguration& configuration){
    Voltage v = _manager.getVoltage(configuration);
    return configuration.numWorkers*configuration.frequency*v*v;
}

void PredictorSimple::prepareForPredictions(){
    ;
}

double PredictorSimple::predict(const FarmConfiguration& configuration){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _manager._samples->average().bandwidth *
                   getScalingFactor(configuration);
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(configuration);
        }break;
    }
    return 0.0;
}

Calibrator::Calibrator(ManagerFarm& manager):
        _manager(manager), _state(CALIBRATION_SEEDS),
        _numCalibrationPoints(0), _calibrationStartMs(0),
        _firstPointGenerated(false){
    _minNumPoints = std::max(_manager._primaryPredictor->getMinimumPointsNeeded(),
                             _manager._secondaryPredictor->getMinimumPointsNeeded());
    DEBUG("Minimum number of points required for calibration: " << _minNumPoints);
    //TODO Assicurarsi che il numero totale di configurazioni possibili sia maggiore del numero minimo di punti
}

void Calibrator::refine(){
    // We have to refine only if a new configuration will be generated.
    if(_state == CALIBRATION_FINISHED &&
       !_manager.isContractViolated()){
        return;
    }
    _manager._primaryPredictor->refine();
    _manager._secondaryPredictor->refine();
}

bool Calibrator::highError() const{
    double primaryValue = _manager.getPrimaryValue();
    double secondaryValue = _manager.getSecondaryValue();
    double primaryError = std::abs((primaryValue - _manager._primaryPrediction)/
                                   primaryValue)*100.0;
    double secondaryError;

    // TODO: This check now works because both service time and power are always  > 0
    // In the future we must find another way to indicat that secondary prediction
    // has not been done.
    if(_manager._secondaryPrediction >= 0.0){
        secondaryError = std::abs((secondaryValue - _manager._secondaryPrediction)/
                                   secondaryValue)*100.0;
    }else{
        secondaryError = 0.0;
    }
    DEBUG("Primary prediction: " << _manager._primaryPrediction << " " <<
          "Secondary prediction: " << _manager._secondaryPrediction);
    DEBUG("Primary error: " << primaryError << " " <<
          "Secondary error: " << secondaryError);
    return primaryError > _manager._p.maxPrimaryPredictionError ||
           secondaryError > _manager._p.maxSecondaryPredictionError;
}

FarmConfiguration Calibrator::getNextConfiguration(){
    FarmConfiguration fc;

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
        refine();
    }else{
        _firstPointGenerated = true;
    }

    switch(_state){
        case CALIBRATION_SEEDS:{
            fc = generateConfiguration();
            if(_numCalibrationPoints + 1 >= _minNumPoints){
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: Moving to predict");
            }
        }break;
        case CALIBRATION_TRY_PREDICT:{
            fc = _manager.getNewConfiguration();
            _state = CALIBRATION_EXTRA_POINT;
            DEBUG("[Calibrator]: Moving to extra");
        }break;
        case CALIBRATION_EXTRA_POINT:{
            if(highError()){
                fc = generateConfiguration();
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: High error");
                DEBUG("[Calibrator]: Moving to predict");
            }else if(_manager.isContractViolated()){
                DEBUG("[Calibrator]: Contract violated");
                fc = _manager.getNewConfiguration();
                _state = CALIBRATION_TRY_PREDICT;
                DEBUG("[Calibrator]: Moving to predict");
            }else{
                fc = _manager.getNewConfiguration();
                _state = CALIBRATION_FINISHED;
                _manager._primaryPredictor->clear();
                _manager._secondaryPredictor->clear();

                CalibrationStats cs;
                // We do -1 because we counted the current point and now we
                // discovered it isn't a calibration point.
                cs.numSteps = _numCalibrationPoints - 1;
                cs.duration = (getMillisecondsTime() - _calibrationStartMs);
                _calibrationStats.push_back(cs);
                DEBUG("[Calibrator]: Moving to finished");
                DEBUG("[Calibrator]: Finished in " << _numCalibrationPoints - 1 <<
                      " steps with configuration " << fc);
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
            if(_manager.isContractViolated()){
                reset();
                fc = generateConfiguration();
                _numCalibrationPoints = 1;
                _state = CALIBRATION_SEEDS;
                _calibrationStartMs = getMillisecondsTime();
                DEBUG("[Calibrator]: Moving to seeds");
            }else{
                fc = _manager._currentConfiguration;
            }
        }break;
    }

    if(_state != CALIBRATION_FINISHED){
        ++_numCalibrationPoints;
    }

    return fc;
}

std::vector<CalibrationStats> Calibrator::getCalibrationsStats() const{
    return _calibrationStats;
}

CalibratorLowDiscrepancy::CalibratorLowDiscrepancy(ManagerFarm& manager):
        Calibrator(manager), _manager(manager){
    uint d = _manager.getConfigurationDimension();
    DEBUG("[Calibrator] Will generate low discrepancy points in " << d <<
          " dimensions");
    const gsl_qrng_type* generatorType;
    switch(_manager._p.strategyCalibration){
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
                                     _manager._p.strategyCalibration);
        }break;
    }
    _generator = gsl_qrng_alloc(generatorType, d);
    _normalizedPoint = new double[d];
}

CalibratorLowDiscrepancy::~CalibratorLowDiscrepancy(){
    gsl_qrng_free(_generator);
    delete[] _normalizedPoint;
}

FarmConfiguration CalibratorLowDiscrepancy::generateConfiguration() const{
    FarmConfiguration r;
    size_t nextPointId = 0;
    gsl_qrng_get(_generator, _normalizedPoint);

    if(_manager.reconfigureWorkers()){
        r.numWorkers =std::ceil(_normalizedPoint[nextPointId++]*
                                (double)_manager._maxNumWorkers);
        if(!r.numWorkers){r.numWorkers = 1;}
    }

    if(_manager.reconfigureFrequency()){
        size_t frequencyId = std::ceil(_normalizedPoint[nextPointId++]*
                                       (double)_manager._availableFrequencies.size());
        if(frequencyId == _manager._availableFrequencies.size()){--frequencyId;}
        r.frequency = _manager._availableFrequencies.at(frequencyId);
    }
    return r;
}

void CalibratorLowDiscrepancy::reset(){
    gsl_qrng_init(_generator);
}

}

#endif //PREDICTORS_IMPL_HPP_
