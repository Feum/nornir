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

#include <gsl/gsl_qrng.h>

using namespace mlpack::regression;

namespace adpff{

class AdaptivityManagerFarm;
class FarmConfiguration;

/**
 * Represents a sample to be used in the regression.
 */
class RegressionData{
public:
    /**
     * Initializes the regression data.
     * @param configuration The configuration to insert as data.
     */
    virtual void init(const FarmConfiguration& configuration) = 0;

    /**
     * Transforms this sample into an Armadillo's row.
     */
    virtual void toArmaRow(size_t rowId, arma::mat& matrix) const = 0;

    /**
     * Gets the number of predictors.
     * @return The number of predictors.
     */
    virtual uint getNumPredictors() const = 0;

    virtual ~RegressionData(){;}
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of bandwidth.
 */
class RegressionDataServiceTime: public RegressionData{
private:
    const AdaptivityManagerFarm& _manager;
    double _physicalCores;
    double _invScalFactorPhysical;
    double _workers;
    double _invScalFactorWorkers;
    double _frequency;

    uint _numPredictors;
public:
    void init(const FarmConfiguration& configuration);

    RegressionDataServiceTime(const AdaptivityManagerFarm& manager);

    RegressionDataServiceTime(const AdaptivityManagerFarm& manager, const FarmConfiguration& configuration);

    uint getNumPredictors() const;

    void toArmaRow(size_t rowId, arma::mat& matrix) const;
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of power.
 */
class RegressionDataPower: public RegressionData{
private:
    const AdaptivityManagerFarm& _manager;
    double _dynamicPowerModel;
    double _voltagePerUsedSockets;
    double _voltagePerUnusedSockets;
    unsigned int _additionalContextes;

    uint _numPredictors;
public:
    void init(const FarmConfiguration& configuration);

    RegressionDataPower(const AdaptivityManagerFarm& manager);

    RegressionDataPower(const AdaptivityManagerFarm& manager, const FarmConfiguration& configuration);

    uint getNumPredictors() const;

    void toArmaRow(size_t rowId, arma::mat& matrix) const;
};


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
     * Gets the number of minimum points needed.
     * Let this number be x. The user need to call
     * refine() method at least x times before
     * starting doing predictions.
     */
    virtual uint getMinimumPointsNeeded(){return 0;}

    /**
     * Clears the predictor removing all the collected
     * data
     */
    virtual void clear(){;}

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

/*
 * Represents a linear regression predictor.
 */
class PredictorLinearRegression: public Predictor{
private:
    PredictorType _type;
    const AdaptivityManagerFarm& _manager;
    LinearRegression _lr;
    std::vector<RegressionData*> _data;
    std::vector<double> _responses; ///< The responses to be used in the regression.
    RegressionData* _predictionInput; ///< Input to be used for predicting a value.
public:
    PredictorLinearRegression(PredictorType type, const AdaptivityManagerFarm& manager);

    ~PredictorLinearRegression();

    void clear();

    uint getMinimumPointsNeeded();

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
    const AdaptivityManagerFarm& _manager;

    double getScalingFactor(const FarmConfiguration& configuration);

    double getPowerPrediction(const FarmConfiguration& configuration);
public:
    PredictorSimple(PredictorType type, const AdaptivityManagerFarm& manager);

    void prepareForPredictions();

    double predict(const FarmConfiguration& configuration);
};

/**
 * State of calibration.
 * Used to track the process.
 */
typedef enum{
    CALIBRATION_SEEDS = 0,
    CALIBRATION_TRY_PREDICT,
    CALIBRATION_EXTRA_POINT,
    CALIBRATION_FINISHED
}CalibrationState;

/**
 * Used to obtain calibration points.
 */
class Calibrator{
public:
    virtual ~Calibrator(){;}

    virtual FarmConfiguration getNextConfiguration() = 0;
};

/**
 * It chooses the calibration points using low discrepancy
 * sampling techniques.
 */
class CalibratorLowDiscrepancy: public Calibrator{
private:
    AdaptivityManagerFarm& _manager;
    CalibrationState _state;
    gsl_qrng* _generator;
    double* _normalizedPoint;
    size_t _minNumPoints;
    size_t _numGeneratedPoints;

    /**
     * Returns true if the error of the model is too high, false otherwise.
     */
    bool highError();

    FarmConfiguration generateConfiguration() const;
public:
    CalibratorLowDiscrepancy(AdaptivityManagerFarm& manager);
    ~CalibratorLowDiscrepancy();

    FarmConfiguration getNextConfiguration();
};

}


#endif /* PREDICTORS_HPP_ */
