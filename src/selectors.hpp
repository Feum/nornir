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
    const FarmConfiguration& _configuration;
    const Smoother<MonitoredSample>* _samples;
    size_t _numCalibrationPoints;

    /**
     * Checks if a specific primary value respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasiblePrimaryValue(double value, bool conservative) const;

    /**
     * Checks if the contract is violated.
     * @param primaryValue The primary value.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated(double primaryValue) const;

    /**
     * Starts the recording of calibration stats.
     * @param totalTasks The total number of tasks processed up to now.
     */
    void startCalibration(uint64_t totalTasks);

    /**
     * Checks if the application phase changed.
     * @param primaryValue The primary value.
     * @param secondaryValue The secondary value.
     * @return true if the phase changed, false otherwise.
     */
    bool phaseChanged(double primaryValue, double secondaryValue) const;
public:
    Selector(const Parameters& p,
             const FarmConfiguration& configuration,
             const Smoother<MonitoredSample>* samples);

    virtual ~Selector(){;}

    /**
     * Returns the next values to be set for the knobs.
     * @param primaryValue The primary value.
     * @param secondaryValue The secondary value.
     * @param totalTasks The total processed tasks.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues(double primaryValue,
                                           double secondaryValue,
                                           u_int64_t totalTasks) = 0;

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
};

/**
 * Always returns the current configuration.
 */
class SelectorFixed: public Selector{
public:
    SelectorFixed(const Parameters& p,
             const FarmConfiguration& configuration,
             const Smoother<MonitoredSample>* samples);

    ~SelectorFixed();

    KnobsValues getNextKnobsValues(double primaryValue,
                                   double secondaryValue,
                                   u_int64_t totalTasks);
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
     * @param primaryValue The primary value.
     * @return The best relative knobs values.
     */
    KnobsValues getBestKnobsValues(double primaryValue);

    /**
     * Refines the models with current data.
     */
    bool refine();

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
     * @param primaryValue The primary value.
     * @param secondaryValue The secondary value.
     * @return True if the predictions were accurate, false otherwise.
     */
    bool isAccurate(double primaryValue, double secondaryValue);
public:
    SelectorPredictive(const Parameters& p,
                       const FarmConfiguration& configuration,
                       const Smoother<MonitoredSample>* samples,
                       std::unique_ptr<Predictor> bandwidthPredictor,
                       std::unique_ptr<Predictor> powerPredictor);

    virtual ~SelectorPredictive();

    /**
     * Returns the next values to be set for the knobs.
     * @param primaryValue The primary value.
     * @param secondaryValue The secondary value.
     * @param totalTasks The total processed tasks.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues(double primaryValue,
                                           double secondaryValue,
                                           u_int64_t totalTasks) = 0;


};

/**
 * Selector described in PDP2015 paper.
 */
class SelectorAnalytical: public SelectorPredictive{
public:
    SelectorAnalytical(const Parameters& p,
                   const FarmConfiguration& configuration,
                   const Smoother<MonitoredSample>* samples);

    KnobsValues getNextKnobsValues(double primaryValue,
                                   double secondaryValue,
                                   u_int64_t totalTasks);
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
    KnobsValues _previousConfiguration;
public:
    SelectorLearner(const Parameters& p,
                       const FarmConfiguration& configuration,
                       const Smoother<MonitoredSample>* samples);

    ~SelectorLearner();

    KnobsValues getNextKnobsValues(double primaryValue,
                                   double secondaryValue,
                                   u_int64_t totalTasks);

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

    mammut::cpufreq::Frequency findNearestFrequency(mammut::cpufreq::Frequency f) const;
    void goRight();
    void goLeft();
public:
    SelectorLiMartinez(const Parameters& p,
                         const FarmConfiguration& configuration,
                         const Smoother<MonitoredSample>* samples);
    ~SelectorLiMartinez();

    KnobsValues getNextKnobsValues(double primaryValue,
                                   double secondaryValue,
                                   u_int64_t totalTasks);
};

}



#endif /* NORNIR_SELECTORS_HPP_ */
