/*
 * predictors.hpp
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

/*! \file predictors.hpp
 * \brief Predictors used by adaptive farm.
 *
 */
#ifndef PREDICTORS_HPP_
#define PREDICTORS_HPP_


#include "farm.hpp"

#include <mlpack/core.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
using namespace mlpack::regression;

namespace adpff{


/**
 * Represents a sample to be used in the regression.
 */
class RegressionData{
public:
    /**
     * Initializes the regression data.
     * @param manager The adaptive farm manager.
     * @param configuration The configuration to insert as data.
     */
    virtual void init(const AdaptivityManagerFarm<>& manager, FarmConfiguration configuration) = 0;

    /**
     * Transforms this sample into an Armadillo's row.
     */
    virtual arma::subview_row toArmaRow() const = 0;

    virtual ~RegressionData(){;}
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of bandwidth.
 */
class RegressionDataBandwidth: public RegressionData{
private:
    unsigned int _physicalCoresInverse;
    unsigned int _workersInverse;
    double _frequencyInverse;
    double _scalFactorPhysical;
    double _scalFactorWorkers;
public:
    void init(const AdaptivityManagerFarm<>& manager, FarmConfiguration configuration){
        double usedPhysicalCores = std::min(configuration.numWorkers, manager._numPhysicalCores);
        _physicalCoresInverse = 1.0 / usedPhysicalCores;
        _workersInverse = 1.0 / (double)configuration.numWorkers;
        _frequencyInverse = 1.0 / (double) configuration.frequency;
        _scalFactorPhysical = ( usedPhysicalCores * ((double)configuration.frequency / (double)manager._availableFrequencies.at(0)));
        _scalFactorWorkers = ((double) configuration.numWorkers * ((double)configuration.frequency / (double)manager._availableFrequencies.at(0)));
    }

    RegressionDataBandwidth():_physicalCoresInverse(0), _workersInverse(0),
                              _frequencyInverse(0), _scalFactorPhysical(0), _scalFactorWorkers(0){;}

    RegressionDataBandwidth(const AdaptivityManagerFarm<>& manager, FarmConfiguration configuration){
        init(manager, configuration);
    }

    arma::subview_row toArmaRow() const{
        arma::subview_row row;
        row(0) = _physicalCoresInverse;
        row(1) = _workersInverse;
        row(2) = _frequencyInverse;
        row(3) = _scalFactorPhysical;
        row(4) = _scalFactorWorkers;
        return row;
    }
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of power.
 */
class RegressionDataPower: public RegressionData{
private:
    unsigned int _contextes;
    double _voltagePerUsedSockets;
    double _voltagePerUnusedSockets;
    double _dynamicPowerModel;
public:
    void init(const AdaptivityManagerFarm<>& manager, FarmConfiguration configuration){
        double usedPhysicalCores = std::min(configuration.numWorkers, manager._numPhysicalCores);
        double voltage = manager.getVoltage(configuration);
        _contextes = std::ceil((double)configuration.numWorkers / usedPhysicalCores);
        _voltagePerUsedSockets = voltage * (double) manager._usedCpus.size();
        _voltagePerUnusedSockets = voltage * (double) manager._unusedCpus.size();
        _dynamicPowerModel = (usedPhysicalCores*configuration.frequency*voltage*voltage)/1000000.0;
    }

    RegressionDataPower():_contextes(0), _voltagePerUsedSockets(0), _voltagePerUnusedSockets(0), _dynamicPowerModel(0){;}

    RegressionDataPower(const AdaptivityManagerFarm<>& manager, FarmConfiguration configuration){
        init(manager, configuration);
    }

    arma::subview_row toArmaRow() const{
        arma::subview_row row;
        row(0) = _contextes;
        row(1) = _voltagePerUsedSockets;
        row(2) = _voltagePerUnusedSockets;
        row(3) = _dynamicPowerModel;
        return row;
    }
};

/*
 * Represents a linear regression predictor.
 */
class PredictorLinearRegression: public Predictor{
private:
    PredictorType _type;
    const AdaptivityManagerFarm<>& _manager;
    LinearRegression _lr;
    RegressionData** _data;
    size_t _dataIndex;
    size_t _dataSize;
    Window<double> _responses; ///< The responses to be used in the regression.
    RegressionData* _predictionInput; ///< Input to be used for predicting a value.
public:
    PredictorLinearRegression(PredictorType type, const AdaptivityManagerFarm<>& manager):
        _type(type), _manager(manager), _dataIndex(0), _dataSize(0),
        _responses(_manager._p.numRegressionPoints){
        _data = new RegressionData*[_manager._p.numRegressionPoints];

        switch(_type){
            case PREDICTION_BANDWIDTH:{
                for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
                    _data[i] = new RegressionDataBandwidth();
                }
                _predictionInput = new RegressionDataBandwidth();
            }break;
            case PREDICTION_POWER:{
                for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
                    _data[i] = new RegressionDataPower();
                }
                _predictionInput = new RegressionDataPower();
            }break;
        }
    }

    ~PredictorLinearRegression(){
        for(size_t i = 0; i < _manager._p.numRegressionPoints; i++){
            delete _data[i];
        }
        delete[] _data;
        delete _predictionInput;
    }

    void refine(){
        _data[_dataIndex]->init(_manager, _manager._currentConfiguration);
        switch(_type){
            case PREDICTION_BANDWIDTH:{
                _responses.add(_manager._averageTasks / (double) _manager._p.samplingInterval); //TODO: Probabilmente bisognerebbe usare il Ts
            }break;
            case PREDICTION_POWER:{
                _responses.add((_manager._usedJoules.cores / _manager._p.samplingInterval) +
                               (_manager._unusedJoules.cores /_manager._p.samplingInterval)); //TODO: Probabilmente si può evitare di dividere per il sampling interval
            }break;
        }
        _dataIndex = (_dataIndex + 1) % _manager._p.numRegressionPoints;
        if(_dataSize < _manager._p.numRegressionPoints){
            ++_dataSize;
        }
    }

    void prepareForPredictions(){
        arma::mat dataMl;
        arma::vec responsesMl;

        for(size_t i = 0; i < _dataSize; i++){
            dataMl.row(i) = _data[i]->toArmaRow();
            responsesMl(i) = _responses[i];
        }

        _lr(dataMl, responsesMl);
    }

    double predict(const FarmConfiguration& configuration){
        _predictionInput->init(_manager, configuration);
        arma::mat predictionInputMl;
        arma::vec result;
        predictionInputMl.row(0) = _predictionInput->toArmaRow();

        _lr.Predict(predictionInputMl, result);
        return result.at(0);
    }
};

/*
 * Represents a simple predictor. It works only for application that
 * exhibits good scalability e do not use hyperthreading.
 * For power, the prediction preserves the relative order between
 * configurations but does not give an exact value. Accordingly,
 * it can't be used for power bounded contracts.
 */
class PredictorSimple: public Predictor{
private:
    PredictorType _type;
    const AdaptivityManagerFarm<>& _manager;
    time_t _now;

    double getScalingFactor(const FarmConfiguration& configuration){
        return (double)(configuration.frequency * configuration.numWorkers) /
               (double)(_manager._currentConfiguration.frequency * _manager._currentConfiguration.numWorkers);
    }

    double getPowerPrediction(const FarmConfiguration& configuration){
        cpufreq::Voltage v = _manager.getVoltage(configuration);
        return configuration.numWorkers*configuration.frequency*v*v;
    }
public:
    PredictorSimple(PredictorType type, const AdaptivityManagerFarm<>& manager):
        _type(type), _manager(manager), _now(0){
        ;
    }

    void prepareForPredictions(){
        _now = time(NULL);
    }

    double predict(const FarmConfiguration& configuration){
        switch(_type){
            case PREDICTION_UTILIZATION:{
                return _manager.getMonitoredValue() * (1.0 / getScalingFactor(configuration));
            }break;
            case PREDICTION_BANDWIDTH:{
                return _manager.getMonitoredValue() * getScalingFactor(configuration);
            }break;
            case PREDICTION_POWER:{
                return getPowerPrediction(configuration);
            }break;
        }
    }
};

};


#endif /* PREDICTORS_HPP_ */
