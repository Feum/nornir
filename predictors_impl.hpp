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

#if 1
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#else
#define DEBUG(x)
#endif


namespace adpff{


/**************** PredictorLinearRegression ****************/

void RegressionDataServiceTime::init(const FarmConfiguration& configuration){
    _numPredictors = 0;

    double usedPhysicalCores = std::min(configuration.numWorkers, _manager._numPhysicalCores);
    /*
    _physicalCores = usedPhysicalCores;
    ++_numPredictors;
     */
    _invScalFactorPhysical = 1.0 / (usedPhysicalCores * ((double)configuration.frequency /
                                                         (double)_manager._availableFrequencies.at(0)));
    ++_numPredictors;

    /*
    if(_manager._p.strategyHyperthreading != STRATEGY_HT_NO){
        _workers = (double)configuration.numWorkers;
        ++_numPredictors;

        _invScalFactorWorkers = 1.0 / ((double) configuration.numWorkers * ((double)configuration.frequency /
                                                                            (double)_manager._availableFrequencies.at(0)));
        ++_numPredictors;
    }
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        _frequency = configuration.frequency;
        ++_numPredictors;
    }*/
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
    //matrix(rowId++, columnId) = _physicalCores;
    matrix(rowId++, columnId) = _invScalFactorPhysical;
    /*if(_manager._p.strategyHyperthreading != STRATEGY_HT_NO){
        matrix(rowId++, columnId) = _workers;
        matrix(rowId++, columnId) = _invScalFactorWorkers;
    }
    if(_manager._p.strategyFrequencies == STRATEGY_FREQUENCY_YES){
        matrix(rowId++, columnId) = _frequency;
    }*/
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
        usedCpus = std::max(usedCpus, _manager._numCpus);
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
        _manager(manager), _dynamicPowerModel(0), _voltagePerUsedSockets(0),
        _voltagePerUnusedSockets(0), _additionalContextes(0), _numPredictors(0){;}

RegressionDataPower::RegressionDataPower(const AdaptivityManagerFarm& manager, const FarmConfiguration& configuration):
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

PredictorLinearRegression::PredictorLinearRegression(PredictorType type, const AdaptivityManagerFarm& manager):
                                                    _type(type), _manager(manager), _dataIndex(0), _dataSize(0),
                                                    _responses(_manager._p.numRegressionPoints){
    _data = new RegressionData*[_manager._p.numRegressionPoints];

    switch(_type){
        case PREDICTION_BANDWIDTH:{
            for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
                _data[i] = new RegressionDataServiceTime(_manager);
            }
            _predictionInput = new RegressionDataServiceTime(_manager);
        }break;
        case PREDICTION_POWER:{
            for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
                _data[i] = new RegressionDataPower(_manager);
            }
            _predictionInput = new RegressionDataPower(_manager);
        }break;
    }
}

PredictorLinearRegression::~PredictorLinearRegression(){
    for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
        delete _data[i];
    }
    delete[] _data;
    delete _predictionInput;
}

void PredictorLinearRegression::refine(){
    _data[_dataIndex]->init(_manager._currentConfiguration);
    double response = 0;
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            response = 1.0 / _manager._averageBandwidth;
        }break;
        case PREDICTION_POWER:{
            response = _manager._averageWatts.cores;
        }break;
    }
    _responses.add(response);
    DEBUG("Refining with configuration [" << _manager._currentConfiguration.numWorkers << ", "
                                          << _manager._currentConfiguration.frequency << "]: "
                                          << response);
    _dataIndex = (_dataIndex + 1) % _manager._p.numRegressionPoints;
    if(_dataSize < _manager._p.numRegressionPoints){
        ++_dataSize;
    }
}

void PredictorLinearRegression::prepareForPredictions(){
    if(!_dataSize){
        return;
    }

    // One observation per column.
    arma::mat dataMl(_data[0]->getNumPredictors(), _dataSize);
    arma::vec responsesMl(_dataSize);
    for(size_t i = 0; i < _dataSize; i++){
        _data[i]->toArmaRow(i, dataMl);
        responsesMl(i) = _responses[i];
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
    DEBUG("Prediction at configuration [" << configuration.numWorkers << ", " << configuration.frequency << "]: " << result.at(0));
    return result.at(0);
}
    
/**************** PredictorSimple ****************/

PredictorSimple::PredictorSimple(PredictorType type, const AdaptivityManagerFarm& manager):
    _type(type), _manager(manager), _now(0){
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
    _now = time(NULL);
}

double PredictorSimple::predict(const FarmConfiguration& configuration){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _manager.getMonitoredValue() * getScalingFactor(configuration);
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(configuration);
        }break;
    }
    return 0.0;
}

CalibratorSpread::CalibratorSpread(const AdaptivityManagerFarm& manager):_manager(manager){;}

std::vector<FarmConfiguration> CalibratorSpread::getCalibrationPoints(){
    /*
     * If I need X points, I split the interval into X - 1 parts and I took
     * the intervals' bounds as points.
     *
     * E.g. X = 4:
     *
     * P1   P2   P3   P4
     * |----|----|----|
     *
     */
    size_t numWorkers = _manager._maxNumWorkers ;
    std::vector<cpufreq::Frequency> frequencies = _manager._availableFrequencies;
    size_t numFrequencies = frequencies.size();
    size_t m = _manager._p.numRegressionPoints - 1;

    std::vector<FarmConfiguration> r;

    size_t nextWorkerId = 1;
    size_t nextFrequencyIndex = 0;

    for(size_t i = 0; i < m; i++){
        r.push_back(FarmConfiguration(nextWorkerId, frequencies.at(nextFrequencyIndex)));
        DEBUG("Calibration point: [" << r.back().numWorkers << ", " << r.back().frequency << "]");

        nextWorkerId += ((numWorkers*i+numWorkers)/m - (numWorkers*i)/m);
        nextFrequencyIndex += ((numFrequencies*i+numFrequencies)/m - (numFrequencies*i)/m);
    }

    r.push_back(FarmConfiguration(numWorkers, frequencies.at(numFrequencies - 1)));
    DEBUG("Calibration point: [" << r.back().numWorkers << ", " << r.back().frequency << "]");
    return r;
}

}

#endif //PREDICTORS_IMPL_HPP_
