/*
 * selectors.hpp
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

#ifndef NORNIR_SELECTORS_HPP_
#define NORNIR_SELECTORS_HPP_

#include "predictors.hpp"
#include "explorers.hpp"
#include "configuration.hpp"
#include "knob.hpp"
#include "external/Mammut/mammut/mammut.hpp"

#include <memory>

namespace nornir{

class Selector{
private:
    std::vector<CalibrationStats> _calibrationStats;
    uint _calibrationStartMs;
    uint64_t _calibrationStartTasks;
    mammut::Mammut _localMammut;
    mammut::energy::Counter* _joulesCounter;
    double _totalCalibrationTime;
    bool _calibrating;
protected:
    const Parameters& _p;
    const Configuration& _configuration;
    const Smoother<MonitoredSample>* _samples;
    size_t _numCalibrationPoints;
    KnobsValues _previousConfiguration;
    Smoother<double>* _bandwidthIn;
    bool _forced;
    bool _forcedReturned;
    KnobsValues _forcedConfiguration;
    bool _calibrationCoordination;
    bool _calibrationAllowed;

    /**
     * Checks if a specific primary value respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasiblePrimaryValue(double value, bool conservative) const;

    /**
     * Checks if the contract is violated.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated() const;

    /**
     * Starts the recording of calibration stats.
     * @param totalTasks The total number of tasks processed up to now.
     */
    void startCalibration(uint64_t totalTasks);
public:
    Selector(const Parameters& p,
             const Configuration& configuration,
             const Smoother<MonitoredSample>* samples);

    virtual ~Selector(){;}


    /**
     * Forces the selector on a specific configuration.
     * @param kv The configuration.
     */
    void forceConfiguration(KnobsValues& kv);

    /**
     * Updates the input bandwidth history with the current value.
     * MUST be called before calling getNextKnobsValues(...).
     */
    void updateBandwidthIn();

    /**
     * Returns the next values to be set for the knobs.
     * @param totalTasks The total processed tasks.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues(u_int64_t totalTasks) = 0;

    /**
     * Returns the calibration statistics.
     * @return A vector containing the calibration statistics.
     */
    std::vector<CalibrationStats> getCalibrationsStats() const;

    /**
     * Returns true if the calibrator is in the calibration phase,
     * false otherwise.
     * @return true if the calibrator is in the calibration phase,
     * false otherwise.
     */
    bool isCalibrating() const;


    /**
     * Stops the recording of calibration stats.
     * @param totalTasks The total number of tasks processed up to now.
     */
    void stopCalibration(uint64_t totalTasks);

    /**
     * Returns the total calibration time.
     * @return The total calibration time.
     */
    double getTotalCalibrationTime() const;

    /**
     * Resets the total calibration time.
     */
    void resetTotalCalibrationTime();

    /**
     * If this function is called, the selector needs to coordinate with a
     * centralised manager before performing calibrations.
     */
    void setCalibrationCoordination();

    /**
     * Allows the selector to start calibration.
     */
    void allowCalibration();
};

/**
 * Always returns the current configuration.
 */
class SelectorFixed: public Selector{
public:
    SelectorFixed(const Parameters& p,
             const Configuration& configuration,
             const Smoother<MonitoredSample>* samples);

    ~SelectorFixed();

    KnobsValues getNextKnobsValues(u_int64_t totalTasks);
};

/**
 * A generic selector that uses the predictions of ALL the configurations
 * to find the best one.
 */
class SelectorPredictive: public Selector{
private:
    std::unique_ptr<Predictor> _primaryPredictor;
    std::unique_ptr<Predictor> _secondaryPredictor;
    bool _feasible;
    // Association between REAL values and observed data.
    std::map<KnobsValues, MonitoredSample> _observedValues;
    std::map<KnobsValues, double> _primaryPredictions;
    std::map<KnobsValues, double> _secondaryPredictions;

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
protected:
    double _primaryPrediction;
    double _secondaryPrediction;

    /**
     * Computes the best relative knobs values for the farm.
     * @return The best relative knobs values.
     */
    KnobsValues getBestKnobsValues();

    /**
     * Refines the models with current data.
     */
    void refine();

    /**
     * Updates the predictions for the next configuration.
     * @param next The next configuration.
     */
    void updatePredictions(const KnobsValues& next);

    /**
     * Returns true if the predictors are ready to perform predictions,
     * false otherwise.
     * @return True if the predictors are ready to perform predictions,
     * false otherwise.
     */
    bool predictorsReady() const;

    /**
     * Returns true if predictions have been performed, false otherwise.
     * @return True if predictions have been performed, false otherwise.
     */
    bool predictionsDone() const;

    /**
     * Returns true if the best solution found is a feasible solution,
     * false otherwise.
     */
    bool isBestSolutionFeasible() const;

    /**
     * Clears the predictors.
     */
    void clearPredictors();


    /**
     * Checks the accuracy of the predictions.
     * @return True if the predictions were accurate, false otherwise.
     */
    bool isAccurate();
public:
    SelectorPredictive(const Parameters& p,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>* samples,
                       std::unique_ptr<Predictor> bandwidthPredictor,
                       std::unique_ptr<Predictor> powerPredictor);

    virtual ~SelectorPredictive();

    /**
     * Returns the next values to be set for the knobs.
     * @param totalTasks The total processed tasks.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues(u_int64_t totalTasks) = 0;

    /**
     * Return the primary prediction for a given configuration.
     * @param values The knobs values.
     * @return The primary prediction for a given configuration.
     */
    double getPrimaryPrediction(KnobsValues values);

    /**
     * Return the secondary prediction for a given configuration.
     * @param values The knobs values.
     * @return The secondary prediction for a given configuration.
     */
    double getSecondaryPrediction(KnobsValues values);

    /**
     * Returns a map with all the primary predictions.
     * @return A map with all the primary predictions.
     */
    const std::map<KnobsValues, double>& getPrimaryPredictions() const;

    /**
     * Returns a map with all the secondary predictions.
     * @return A map with all the secondary predictions.
     */
    const std::map<KnobsValues, double>& getSecondaryPredictions() const;
};

/**
 * Selector described in PDP2015 paper.
 */
class SelectorAnalytical: public SelectorPredictive{
private:
    uint _violations;
public:
    SelectorAnalytical(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    KnobsValues getNextKnobsValues(u_int64_t totalTasks);
};

/**
 * A generic online learner selector.
 */
class SelectorLearner: public SelectorPredictive{
private:
    Explorer* _explorer;
    bool _firstPointGenerated;
    uint _contractViolations;
    uint _accuracyViolations;
    std::unique_ptr<Predictor> getPredictor(PredictorType type,
                                            const Parameters& p,
                                            const Configuration& configuration,
                                            const Smoother<MonitoredSample>* samples) const;
public:
    SelectorLearner(const Parameters& p,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>* samples);

    ~SelectorLearner();

    KnobsValues getNextKnobsValues(u_int64_t totalTasks);

    /**
     * Checks if the application phase changed.
     * @return true if the phase changed, false otherwise.
     */
    bool phaseChanged() const;
};

/**
 * Explores a fixed number of configurations before starting
 * making predictions.
 * Configuration is never changed once it has been found.
 */
class SelectorFixedExploration: public SelectorPredictive{
private:
    std::vector<KnobsValues> _confToExplore;
public:
    SelectorFixedExploration(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples,
                   std::unique_ptr<Predictor> bandwidthPredictor,
                   std::unique_ptr<Predictor> powerPredictor,
                   size_t numSamples);

    ~SelectorFixedExploration();


    KnobsValues getNextKnobsValues(u_int64_t totalTasks);
};

/**
 * Applies the algorithm described in:
 * "A Probabilistic Graphical Model-based Approach for Minimizing
 * Energy Under Performance Constraints" - Mishra, Nikita and Zhang, Huazhe
 * and Lafferty, John D. and Hoffmann, Henry
 */
class SelectorMishra: public SelectorFixedExploration{
public:
    SelectorMishra(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorMishra();
};

/**
 * This selector does a full exploration of the search space.
 */
class SelectorFullSearch: public SelectorFixedExploration{
public:
    SelectorFullSearch(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorFullSearch();
};

/**
 * A selector that implements the algorithm described in:
 * "Dynamic Power-Performance Adaptation of Parallel Computation
 * on Chip Multiprocessors" - Jian Li and Jose F. Martınez
 */
class SelectorLiMartinez: public Selector{
private:
    bool _firstPointGenerated;
    uint _low1, _mid1, _high1;
    uint _low2, _mid2, _high2;
    uint _midId;
    std::vector<mammut::cpufreq::Frequency> _availableFrequencies;
    double _currentWatts;
    double _optimalWatts;
    mammut::cpufreq::Frequency _optimalFrequency;
    uint _optimalWorkers;
    double _currentBw, _leftBw, _rightBw;
    KnobsValues _optimalKv;
    bool _improved;
    std::vector<double> _allowedCores;

    mammut::cpufreq::Frequency findNearestFrequency(mammut::cpufreq::Frequency f) const;
    void goRight();
    void goLeft();
public:
    SelectorLiMartinez(const Parameters& p,
                         const Configuration& configuration,
                         const Smoother<MonitoredSample>* samples);
    ~SelectorLiMartinez();

    KnobsValues getNextKnobsValues(u_int64_t totalTasks);
};

}



#endif /* NORNIR_SELECTORS_HPP_ */
