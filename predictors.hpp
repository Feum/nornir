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

#include "utils.hpp"

#include <ff/farm.hpp>

#include <mlpack/core.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>

using namespace mlpack::regression;

namespace adpff{

template<typename lb_t, typename gt_t>
class AdaptivityManagerFarm;

class FarmConfiguration;


/**
 * Type of predictor.
 */
typedef enum PredictorType{
    PREDICTION_BANDWIDTH = 0,
    PREDICTION_POWER
}PredictorType;

/*!
 * Represents a generic predictor.
 */
class Predictor{
public:
    virtual ~Predictor(){;}

    /**
     * If possible, refines the model with the information obtained on the current
     * configuration.
     */
    virtual void refine(){;}

    /**
     * Prepare the predictor to accept a set of prediction requests.
     */
    virtual void prepareForPredictions() = 0;

    /**
     * Predicts the value at a specific configuration.
     * @param configuration The configuration.
     * @return The predicted value at a specific configuration.
     */
    virtual double predict(const FarmConfiguration& configuration) = 0;
};

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
    virtual void init(const AdaptivityManagerFarm<>& manager, const FarmConfiguration& configuration) = 0;

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
    void init(const AdaptivityManagerFarm<>& manager, const FarmConfiguration& configuration);

    RegressionDataBandwidth();

    RegressionDataBandwidth(const AdaptivityManagerFarm<>& manager, const FarmConfiguration& configuration);

    arma::subview_row toArmaRow() const;
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
    void init(const AdaptivityManagerFarm<>& manager, const FarmConfiguration& configuration);

    RegressionDataPower();

    RegressionDataPower(const AdaptivityManagerFarm<>& manager, const FarmConfiguration& configuration);

    arma::subview_row toArmaRow() const;
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
    PredictorLinearRegression(PredictorType type, const AdaptivityManagerFarm<>& manager);

    ~PredictorLinearRegression();

    void refine();

    void prepareForPredictions();

    double predict(const FarmConfiguration& configuration);
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

    double getScalingFactor(const FarmConfiguration& configuration);

    double getPowerPrediction(const FarmConfiguration& configuration);
public:
    PredictorSimple(PredictorType type, const AdaptivityManagerFarm<>& manager);

    void prepareForPredictions();

    double predict(const FarmConfiguration& configuration);
};

};


#endif /* PREDICTORS_HPP_ */
