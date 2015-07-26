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
#include "farm.hpp"

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

void RegressionDataServiceTime::init(const FarmConfiguration& configuration){
    _numPredictors = 0;
    double usedPhysicalCores = std::min(configuration.numWorkers, _manager._numPhysicalCores);

    _physicalCores = usedPhysicalCores;
    ++_numPredictors;

    _invScalFactorPhysical = 1.0 / (usedPhysicalCores * ((double)configuration.frequency /
                                                         (double)_manager._availableFrequencies.at(0)));
    ++_numPredictors;


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

RegressionDataServiceTime::RegressionDataServiceTime(const AdaptivityManagerFarm& manager):
                                                   _manager(manager),
                                                   _physicalCores(0), _invScalFactorPhysical(0),
                                                   _workers(0), _invScalFactorWorkers(0),
                                                   _frequency(0),
                                                   _numPredictors(0){;}

RegressionDataServiceTime::RegressionDataServiceTime(const AdaptivityManagerFarm& manager, const FarmConfiguration& configuration):
            _manager(manager){
    init(configuration);
}

uint RegressionDataServiceTime::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataServiceTime::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    //TODO: !!!
    matrix(rowId++, columnId) = _physicalCores;
    matrix(rowId++, columnId) = _invScalFactorPhysical;
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

RegressionDataPower::RegressionDataPower(const AdaptivityManagerFarm& manager):
        _manager(manager), _dynamicPowerModel(0),
        _voltagePerUsedSockets(0), _voltagePerUnusedSockets(0),
        _additionalContextes(0), _numPredictors(0){;}

RegressionDataPower::RegressionDataPower(const AdaptivityManagerFarm& manager,
                                         const FarmConfiguration& configuration):
            _manager(manager){
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
                                                     const AdaptivityManagerFarm& manager):
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
    for(size_t i = 0; i < _data.size(); i++){
        delete _data.at(i);
    }
    _data.clear();
    _responses.clear();
}

uint PredictorLinearRegression::getMinimumPointsNeeded(){
    //TODO: Brutto a vedersi, ripulire
    _predictionInput->init(_manager._currentConfiguration);
    return std::max<uint>(_predictionInput->getNumPredictors(), 2);
}

void PredictorLinearRegression::refine(){
    RegressionData* rd = NULL;
    double response = 0;

    switch(_type){
        case PREDICTION_BANDWIDTH:{
            rd = new RegressionDataServiceTime(_manager, _manager._currentConfiguration);
            response = 1.0 / _manager._averageBandwidth;
        }break;
        case PREDICTION_POWER:{
            rd = new RegressionDataPower(_manager, _manager._currentConfiguration);
            response = _manager._averageWatts.cores;
        }break;
    }

    _data.push_back(rd);
    _responses.push_back(response);
    DEBUG("Refining with configuration [" << _manager._currentConfiguration.numWorkers << ", "
                                          << _manager._currentConfiguration.frequency << "]: "
                                          << response);
}

void PredictorLinearRegression::prepareForPredictions(){
    if(!_data.size()){
        return;
    }

    // One observation per column.
    arma::mat dataMl(_data[0]->getNumPredictors(), _data.size());
    arma::vec responsesMl(_data.size());
    for(size_t i = 0; i < _data.size(); i++){
        _data.at(i)->toArmaRow(i, dataMl);
        responsesMl(i) = _responses.at(i);
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
                                 const AdaptivityManagerFarm& manager):
    _type(type), _manager(manager){
    ;
}

double PredictorSimple::getScalingFactor(const FarmConfiguration& configuration){
    return (double)(configuration.frequency * configuration.numWorkers) /
           (double)(_manager._currentConfiguration.frequency * _manager._currentConfiguration.numWorkers);
}

double PredictorSimple::getPowerPrediction(const FarmConfiguration& configuration){
    cpufreq::Voltage v = _manager.getVoltage(configuration);
    return configuration.numWorkers*configuration.frequency*v*v;
}

void PredictorSimple::prepareForPredictions(){
    ;
}

double PredictorSimple::predict(const FarmConfiguration& configuration){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _manager._averageBandwidth * getScalingFactor(configuration);
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(configuration);
        }break;
    }
    return 0.0;
}

CalibratorLowDiscrepancy::CalibratorLowDiscrepancy(AdaptivityManagerFarm& manager):
        _manager(manager), _state(CALIBRATION_SEEDS), _numGeneratedPoints(0){
    uint d = _manager.getConfigurationDimension();
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
    _minNumPoints = std::max(_manager._primaryPredictor->getMinimumPointsNeeded(),
                             _manager._secondaryPredictor->getMinimumPointsNeeded());
}

CalibratorLowDiscrepancy::~CalibratorLowDiscrepancy(){
    gsl_qrng_free(_generator);
    delete[] _normalizedPoint;
}

bool CalibratorLowDiscrepancy::highError(){
    double primaryValue = _manager.getPrimaryValue();
    double secondaryValue = _manager.getSecondaryValue();
    double primaryError = std::abs((primaryValue - _manager._primaryPrediction)/
                                   primaryValue)*100.0;
    double secondaryError = std::abs((secondaryValue - _manager._secondaryPrediction)/
                                     secondaryValue)*100.0;
    DEBUG("Primary prediction: " << _manager._primaryPrediction << " " <<
          "Secondary prediction: " << _manager._secondaryPrediction);
    DEBUG("Primary error: " << primaryError << " " <<
          "Secondary error: " << secondaryError);
    return primaryError > _manager._p.maxPredictionError ||
           secondaryError > _manager._p.maxPredictionError;
}

FarmConfiguration CalibratorLowDiscrepancy::generateConfiguration() const{
    FarmConfiguration r;
    size_t nextPointId = 0;
    gsl_qrng_get(_generator, _normalizedPoint);

    if(_manager.reconfigureWorkers()){
        r.numWorkers = _normalizedPoint[nextPointId++]*
                       (double)_manager._maxNumWorkers;
        if(!r.numWorkers){r.numWorkers = 1;}
    }

    if(_manager.reconfigureFrequency()){
        size_t frequencyId = _normalizedPoint[nextPointId++]*
                             (double)_manager._availableFrequencies.size();
        if(frequencyId == _manager._availableFrequencies.size()){--frequencyId;}
        r.frequency = _manager._availableFrequencies.at(frequencyId);
    }
    return r;
}

FarmConfiguration CalibratorLowDiscrepancy::getNextConfiguration(){
    FarmConfiguration fc;

    if(_state != CALIBRATION_FINISHED){
        ++_numGeneratedPoints;
    }

    switch(_state){
        case CALIBRATION_SEEDS:{
            fc = generateConfiguration();
            if(_numGeneratedPoints >= _minNumPoints){
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
                DEBUG("[Calibrator]: Moving to finished");
                DEBUG("[Calibrator]: Finished in " << _numGeneratedPoints - 1 <<
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
                gsl_qrng_init(_generator);
                fc = generateConfiguration();
                _numGeneratedPoints = 1;
                _state = CALIBRATION_SEEDS;                                                                                                                                    
                DEBUG("[Calibrator]: Moving to seeds");
            }else{
                fc = _manager._currentConfiguration;
            }
        }break;
    }
    return fc;
}

}

#endif //PREDICTORS_IMPL_HPP_
