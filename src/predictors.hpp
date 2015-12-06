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

#include "configuration.hpp"
#include "utils.hpp"

#include <mammut/mammut.hpp>

#include <ff/farm.hpp>

#include <mlpack/core.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>

#include <gsl/gsl_qrng.h>

#include <map>

using namespace mlpack::regression;

namespace adpff{

class KnobsValues;

/**
 * Represents a sample to be used in the regression.
 */
class RegressionData{
public:
    RegressionData(const ManagerFarm& manager):
        _manager(manager){
        ;
    }

    /**
     * Initializes the regression data with the current configuration.
     */
    void init();

    /**
     * Initializes the regression data.
     * @param values The values to insert as data.
     */

    virtual void init(const KnobsValues& values) = 0;

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
    mammut::cpufreq::Frequency _minFrequency;
    double _invScalFactorFreq;
    double _invScalFactorFreqAndCores;

    uint _numPredictors;
public:
    void init(const KnobsValues& values);

    RegressionDataServiceTime(const ManagerFarm& manager);

    RegressionDataServiceTime(const ManagerFarm& manager,
                              const KnobsValues& values);

    uint getNumPredictors() const;

    void toArmaRow(size_t rowId, arma::mat& matrix) const;
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of power.
 */
class RegressionDataPower: public RegressionData{
private:
    double _dynamicPowerModel;
    double _voltagePerUsedSockets;
    double _voltagePerUnusedSockets;
    unsigned int _additionalContextes;

    uint _numPredictors;
public:
    void init(const KnobsValues& values);

    RegressionDataPower(const ManagerFarm& manager);

    RegressionDataPower(const ManagerFarm& manager,
                        const KnobsValues& values);

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
     * If possible, refines the model with the information
     * obtained on the current configuration.
     */
    virtual void refine(){;}

    /**
     * Prepare the predictor to accept a set of prediction requests.
     */
    virtual void prepareForPredictions() = 0;

    /**
     * Predicts the value at specific knobs values.
     * @param values The values.
     * @return The predicted value at a specific combination of real knobs values.
     */
    virtual double predict(const KnobsValues& realValues) = 0;
};


typedef struct{
    RegressionData* data;
    double response;
}Observation;

/*
 * Represents a linear regression predictor.
 */
class PredictorLinearRegression: public Predictor{
private:
    PredictorType _type;
    LinearRegression _lr;

    typedef std::map<KnobsValues, Observation> Observations;
    Observations _observations;
    typedef Observations::iterator obs_it;

    // Input to be used for predicting a value.
    RegressionData* _predictionInput;

    double getCurrentResponse() const;
public:
    PredictorLinearRegression(PredictorType type, const ManagerFarm& manager);

    ~PredictorLinearRegression();

    void clear();

    uint getMinimumPointsNeeded();

    void refine();

    void prepareForPredictions();

    double predict(const KnobsValues& configuration);
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

    double getScalingFactor(const KnobsValues& values);

    double getPowerPrediction(const KnobsValues& values);
public:
    PredictorSimple(PredictorType type, const ManagerFarm& manager);

    void prepareForPredictions();

    double predict(const KnobsValues& configuration);
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

typedef struct{
    uint numSteps;
    uint duration;
}CalibrationStats;

/**
 * Used to obtain calibration points.
 */
class Calibrator{
private:
    const Parameters& _p;
    const FarmConfiguration& _configuration;
    CalibrationState _state;
    std::vector<CalibrationStats> _calibrationStats;
    size_t _minNumPoints;
    size_t _numCalibrationPoints;
    uint _calibrationStartMs;
    bool _firstPointGenerated;
    mammut::cpufreq::VoltageTable _voltageTable;

    // The predictor of the primary value.
    Predictor* _primaryPredictor;

    // The predictor of the secondary value.
    Predictor* _secondaryPredictor;

    // The prediction done for the primary value for the chosen configuration.
    double _primaryPrediction;

    // The prediction done for the secondary value for the chosen configuration.
    double _secondaryPrediction;

    bool highError(double primaryValue, double secondaryValue) const;
    void refine(bool isContractViolated);

    /**
     * Returns the voltage at a specific combinations of knobs values.
     * @param values The combination.
     * @return The voltage at that combination.
     */
    double getVoltage(const KnobsValues& values) const;

    /**
     * Checks if x is a best suboptimal monitored value than y.
     * @param x The first monitored value.
     * @param y The second monitored value.
     * @return True if x is a best suboptimal monitored value than y,
     *         false otherwise.
     */
    bool isBestSuboptimalValue(double x, double y) const;

    /**
     * Returns true if x is a best secondary value than y, false otherwise.
     */
    bool isBestSecondaryValue(double x, double y) const;

    /**
     * Checks if a specific primary value respects the required contract.
     * @param value The value to be checked.
     * @param tolerance The percentage of tolerance allowed for the check
     */
    bool isFeasiblePrimaryValue(double value, double tolerance = 0) const;

    /**
     * Computes the best relative knobs values for the farm.
     * @return The best relative knobs values.
     */
    KnobsValues getBestKnobsValues(double primaryValue, double secondaryValue);
protected:
    /**
     *  Override this method to provide custom ways to generate
     *  relative knobs values for calibration.
     *  @return The relative knobs values.
     **/
    virtual KnobsValues generateRelativeKnobsValues() const = 0;
    virtual void reset(){;}
public:
    Calibrator(const Parameters& p, const FarmConfiguration& configuration);

    virtual ~Calibrator(){;}

    /**
     * Returns the next values to be set for the knobs.
     * @param isContractViolated True if the contract has been violated,
     *        false otherwise.
     * @param primaryValue The primary value.
     * @param secondaryValue The secondary value.
     *
     * @return The next values to be set for the knobs.
     */
    KnobsValues getNextKnobsValues(bool isContractViolated,
                                   double primaryValue, double secondaryValue);

    std::vector<CalibrationStats> getCalibrationsStats() const;
};

/**
 * It chooses the calibration points using low discrepancy
 * sampling techniques.
 */
class CalibratorLowDiscrepancy: public Calibrator{
private:
    const Parameters& _p;
    const FarmConfiguration& _configuration;
    gsl_qrng* _generator;
    double* _normalizedPoint;
protected:
    KnobsValues generateRelativeKnobsValues() const;
    void reset();
public:
    CalibratorLowDiscrepancy(const Parameters& p,
                             const FarmConfiguration& configuration);
    ~CalibratorLowDiscrepancy();
};

}


#endif /* PREDICTORS_HPP_ */
