/*
 * predictors.hpp
 *
 * Created on: 10/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

/*! \file predictors.hpp
 * \brief Predictors used by adaptive farm.
 *
 */
#ifndef NORNIR_PREDICTORS_HPP_
#define NORNIR_PREDICTORS_HPP_

#include "configuration.hpp"
#include "utils.hpp"

#include "external/Mammut/mammut/mammut.hpp"

#include <gsl/gsl_multifit.h>
#include <mlpack/core.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include "external/leo/leo.h" // Must be included after mlpack

#include <map>

namespace nornir{

class KnobsValues;

/**
 * Represents a sample to be used in the regression.
 */
class RegressionData{
protected:
    const Parameters& _p;
    const Configuration& _configuration;
    const Smoother<MonitoredSample>* _samples;
    mammut::topology::Topology* _topology;
    uint _cpus;
    uint _phyCores;
    uint _phyCoresPerCpu;
    uint _virtCoresPerPhyCores;

    double getUsedPhysicalCores(double numWorkers);
public:
    RegressionData(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

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
    uint _phyCores;
    mammut::cpufreq::Frequency _minFrequency;
    double _invScalFactorFreq;
    double _invScalFactorFreqAndCores;

    uint _numPredictors;
public:
    void init(const KnobsValues& values);

    RegressionDataServiceTime(const Parameters& p,
                              const Configuration& configuration,
                              const Smoother<MonitoredSample>* samples);

    uint getNumPredictors() const;

    void toArmaRow(size_t rowId, arma::mat& matrix) const;
};

/**
 * A row of the data matrix to be used in linear
 * regression for prediction of power.
 */
class RegressionDataPower: public RegressionData{
private:
    StrategyUnusedVirtualCores _strategyUnused;

    double _dynamicPowerModel;
    double _voltagePerUsedSockets;
    double _voltagePerUnusedSockets;
    unsigned int _additionalContextes;
    mammut::cpufreq::Frequency _lowestFrequency;

    uint _numPredictors;
public:
    void init(const KnobsValues& values);

    RegressionDataPower(const Parameters& p,
                        const Configuration& configuration,
                        const Smoother<MonitoredSample>* samples);

    uint getNumPredictors() const;

    void toArmaRow(size_t rowId, arma::mat& matrix) const;

    /**
     * Set the position of the inactive power parameter into the parameters
     * vector into the variable passed as parameter.
     * @param pos The position of the inactive power parameter into the parameters
     * vector.
     * @return True if there inactive power is considered by the model, false
     * otherwise (e.g. if only one CPU is present on the machine).
     */
    bool getInactivePowerPosition(size_t& pos) const;
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
protected:
    PredictorType _type;
    const Parameters& _p;
    const Configuration& _configuration;
    const Smoother<MonitoredSample>* _samples;
    double _modelError;

    double getMaximumBandwidth() const;
    double getCurrentPower() const;
public:
    Predictor(PredictorType type,
              const Parameters& p,
              const Configuration& configuration,
              const Smoother<MonitoredSample>* samples);

    virtual ~Predictor();

    /**
     * Returns true if the predictor is ready to make prediction.
     * If false is returned, we need to refine the model with
     * more points.
     * @return true if the predictor is ready to make prediction,
     *         false otherwise.
     */
    virtual bool readyForPredictions() = 0;

    /**
     * Clears the predictor removing all the collected
     * data
     */
    virtual void clear() = 0;

    /**
     * If possible, refines the model with the information
     * obtained on the current configuration.
     */
    virtual void refine() = 0;

    /**
     * Prepare the predictor to accept a set of prediction requests.
     * ATTENTION: If it is already ready to perform predictions, nothing should
     *            be done.
     */
    virtual void prepareForPredictions() = 0;

    /**
     * Predicts the value at specific knobs values.
     * @param values The values.
     * @return The predicted value at a specific combination of real knobs values.
     */
    virtual double predict(const KnobsValues& realValues) = 0;

    /**
     * Returns the model error, i.e. the error between the observations and the
     * predictions.
     * @return The model error.
     */
    double getModelError(){return _modelError;}
};


typedef struct{
    RegressionData* data;
    double response;
}Observation;

/*
 * A linear regression predictor for <Cores, Frequency> configurations.
 */
class PredictorLinearRegression: public Predictor{
private:
    mlpack::regression::LinearRegression _lr;

    typedef std::map<KnobsValues, Observation> Observations;
    Observations _observations;
    typedef Observations::iterator obs_it;

    // Aging vector, it contains the last regressionAging KnobsValues
    std::vector<KnobsValues> _agingVector;
    size_t _currentAgingId;

    // Input to be used for predicting a value.
    RegressionData* _predictionInput;

    bool _preparationNeeded;

    // Identifier of rows that have been removed to let the dataMl matrix
    // invertible.
    std::vector<size_t> _removedRows;

    uint _otherApplicationsCores;

    double getCurrentResponse() const;
public:
    PredictorLinearRegression(PredictorType type,
                              const Parameters& p,
                              const Configuration& configuration,
                              const Smoother<MonitoredSample>* samples);

    ~PredictorLinearRegression();

    void clear();

    bool readyForPredictions();

    void refine();

    void prepareForPredictions();

    double predict(const KnobsValues& configuration);

    double getInactivePowerParameter() const;
};

/**
 * A regression predictor for <Cores, Frequency, Mapping> configurations.
 * It works by using one predictor for each mapping.
 */
template <class P>
class PredictorRegressionMapping: public Predictor{
private:
    Predictor* _predictors[MAPPING_TYPE_NUM];
public:
    PredictorRegressionMapping(PredictorType type,
                              const Parameters& p,
                              const Configuration& configuration,
                              const Smoother<MonitoredSample>* samples):
                            	  Predictor(type, p, configuration, samples){
    	for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
    		_predictors[i] = new P(type, p, configuration, samples);
    	}
	}

    ~PredictorRegressionMapping(){
        for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
            delete _predictors[i];
        }
    }

    void clear(){
        for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
            _predictors[i]->clear();
        }
    }

    bool readyForPredictions(){
        for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
            if(!_predictors[i]->readyForPredictions()){
                return false;
            }
        }
        return true;
    }

    void refine(){
    	_predictors[(MappingType) _configuration.getRealValue(KNOB_TYPE_MAPPING)]->refine();
    }

    void prepareForPredictions(){
        for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
            _predictors[i]->prepareForPredictions();
        }
    }

    double predict(const KnobsValues& configuration){
        double realValue = 0.0;
        if(configuration.areRelative()){
            _configuration.getKnob(KNOB_TYPE_MAPPING)->getRealFromRelative(configuration[KNOB_TYPE_MAPPING], realValue);
        }else{
            realValue = configuration[KNOB_TYPE_MAPPING];
        }
        return _predictors[(MappingType) realValue]->predict(configuration);
    }
};

/**
 * A predictor that uses Universal Scalability Law to predict application
 * performance.
 */
#define POLYNOMIAL_DEGREE_USL 3 // Second degree polynomial

class PredictorUsl: public Predictor{
private:
    std::vector<double> _xs;
    std::vector<double> _ys;
    gsl_multifit_linear_workspace *_ws;
    gsl_matrix *_cov, *_x;
    gsl_vector *_y, *_c;
    double _chisq;
    std::vector<double> _coefficients;
    bool _preparationNeeded;
    double _maxFreqBw;
    double _minFreqBw;
    double _minFreqCoresBw;
    double _minFrequency;
    double _maxFrequency;
    double _maxCores;
public:
    PredictorUsl(PredictorType type,
                 const Parameters& p,
                 const Configuration& configuration,
                 const Smoother<MonitoredSample>* samples);

    ~PredictorUsl();

    void clear();

    bool readyForPredictions();

    void refine();

    void prepareForPredictions();

    double predict(const KnobsValues& configuration);
};

/*
 * Represents a simple predictor. It works only for application that
 * exhibits good scalability and do not use hyperthreading.
 * For power, the prediction preserves the relative order between
 * configurations but does not give an exact value. Accordingly,
 * it can't be used for power bounded contracts.
 */
class PredictorAnalytical: public Predictor{
private:
    uint _phyCores;
    uint _phyCoresPerCpu;
    uint _cpus;

    double getScalingFactor(const KnobsValues& values);
    double getPowerPrediction(const KnobsValues& values);
public:
    PredictorAnalytical(PredictorType type,
                    const Parameters& p,
                    const Configuration& configuration,
                    const Smoother<MonitoredSample>* samples);

    bool readyForPredictions();

    void prepareForPredictions();

    double predict(const KnobsValues& values);

    void refine();

    void clear();
};

/**
 * Applies the algorithm described in:
 * "A Probabilistic Graphical Model-based Approach for Minimizing
 * Energy Under Performance Constraints" - Mishra, Nikita and Zhang, Huazhe
 * and Lafferty, John D. and Hoffmann, Henry
 */
class PredictorMishra: public Predictor{
private:
    std::map<KnobsValues, size_t> _confIndexes;
    arma::vec _values;
    arma::vec _predictions;
    bool _preparationNeeded;
    size_t _appId;
public:
    PredictorMishra(PredictorType type,
              const Parameters& p,
              const Configuration& configuration,
              const Smoother<MonitoredSample>* samples);

    ~PredictorMishra();

    bool readyForPredictions();

    void clear();

    void refine();

    void prepareForPredictions();

    double predict(const KnobsValues& realValues);
};

/**
 * Applies a full search strategy in order to find
 * the best configuration.
 */
class PredictorFullSearch: public Predictor{
private:
    const std::vector<KnobsValues>& _allConfigurations;
    std::map<KnobsValues, double> _values;
public:
    PredictorFullSearch(PredictorType type,
              const Parameters& p,
              const Configuration& configuration,
              const Smoother<MonitoredSample>* samples);

    ~PredictorFullSearch();

    bool readyForPredictions();

    void clear();

    void refine();

    void prepareForPredictions();

    double predict(const KnobsValues& realValues);
};

}


#endif /* NORNIR_PREDICTORS_HPP_ */
