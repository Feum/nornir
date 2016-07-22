/*
 * selectors.cpp
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

#include "selectors.hpp"
#include "utils.hpp"
#include <cfloat>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_SELECTORS
#define DEBUG(x) do { std::cerr << "[Selectors] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

#define NOT_VALID DBL_MIN

using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::utils;
using namespace mammut::topology;
using namespace std;

namespace nornir{

Selector::Selector(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples):
        _calibrationStartMs(0),
        _calibrationStartTasks(0),
        _totalCalibrationTime(0),
        _calibrating(false),
        _p(p),
        _configuration(configuration),
        _samples(samples),
        _numCalibrationPoints(0){
    _joulesCounter = _localMammut.getInstanceEnergy()->getCounter();
    //TODO Fare meglio con mammut
    //TODO Assicurarsi che il numero totale di configurazioni possibili sia maggiore del numero minimo di punti
#if 0
    // Input bandwidth smoother                                                                                                                                                                                      
    if(p.strategySmoothing == STRATEGY_SMOOTHING_EXPONENTIAL){
        _bandwidthIn = new MovingAverageExponential<double>(p.smoothingFactor);
    }else{
        _bandwidthIn = new MovingAverageSimple<double>(p.smoothingFactor);
    }
#endif
    _bandwidthIn = new MovingAverageSimple<double>(3);
}


bool Selector::isFeasiblePrimaryValue(double value, bool conservative) const{
    double conservativeOffset = 0;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            if(conservative && _p.conservativeValue){
                conservativeOffset = (_p.overloadThresholdFarm - _p.underloadThresholdFarm) *
                                     (_p.conservativeValue / 100.0) / 2.0;
            }
            return value > _p.underloadThresholdFarm + conservativeOffset &&
                   value < _p.overloadThresholdFarm - conservativeOffset;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            if(conservative && _p.conservativeValue){
                conservativeOffset = _p.requiredBandwidth * (_p.conservativeValue / 100.0);
            }
            return value > _p.requiredBandwidth + conservativeOffset;
        }break;
        case CONTRACT_POWER_BUDGET:{
            if(conservative && _p.conservativeValue){
                conservativeOffset = _p.powerBudget * (_p.conservativeValue / 100.0);
            }

            return value < _p.powerBudget - conservativeOffset;
        }break;
        default:{
            return false;
        }break;
    }
    return false;
}

bool Selector::isContractViolated() const{
    double primaryValue = 0.0;
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            primaryValue = _samples->average().utilisation;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            primaryValue = _samples->average().bandwidth;
        }break;
        case CONTRACT_POWER_BUDGET:{
            primaryValue = _samples->average().watts;
        }break;
        default:{
            return 0;
        }break;
    }
    return !isFeasiblePrimaryValue(primaryValue, false);
}

void Selector::startCalibration(uint64_t totalTasks){
    _calibrating = true;
    _numCalibrationPoints = 0;
    _calibrationStartMs = getMillisecondsTime();
    _calibrationStartTasks = totalTasks;
    if(_joulesCounter){
        _joulesCounter->reset();
    }
}

void Selector::stopCalibration(uint64_t totalTasks){
    _calibrating = false;
    if(_numCalibrationPoints){
        CalibrationStats cs;
        cs.numSteps = _numCalibrationPoints;
        cs.duration = (getMillisecondsTime() - _calibrationStartMs);
        _totalCalibrationTime += cs.duration;
        cs.numTasks = totalTasks - _calibrationStartTasks;
        if(_joulesCounter){
            cs.joules = ((CounterCpus*) _joulesCounter)->getJoulesCoresAll();
        }
        _calibrationStats.push_back(cs);
        _numCalibrationPoints = 0;
    }
}

double Selector::getTotalCalibrationTime() const{
    return _totalCalibrationTime;
}

void Selector::resetTotalCalibrationTime(){
    _totalCalibrationTime = 0;
}


std::vector<CalibrationStats> Selector::getCalibrationsStats() const{
    return _calibrationStats;
}

bool Selector::isCalibrating() const{
    return _calibrating;
}

void Selector::updateBandwidthIn(){
    if(_samples->average().utilisation < MAX_RHO){
        if(_bandwidthIn->average() == numeric_limits<double>::max()){
            _bandwidthIn->reset();
        }
        _bandwidthIn->add(_samples->getLastSample().bandwidth);
    }else if(!_bandwidthIn->size()){
        _bandwidthIn->add(numeric_limits<double>::max());
    }
}


SelectorFixed::SelectorFixed(const Parameters& p,
         const Configuration& configuration,
         const Smoother<MonitoredSample>* samples):
                 Selector(p, configuration, samples){
    ;
}

SelectorFixed::~SelectorFixed(){;}

KnobsValues SelectorFixed::getNextKnobsValues(u_int64_t totalTasks){
    _previousConfiguration = _configuration.getRealValues();
    return _configuration.getRealValues();
}

SelectorPredictive::SelectorPredictive(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples,
                   std::unique_ptr<Predictor> bandwidthPredictor,
                   std::unique_ptr<Predictor> powerPredictor):
                       Selector(p, configuration, samples),
                       _feasible(true),
                       _primaryPrediction(NOT_VALID),
                       _secondaryPrediction(NOT_VALID){
    /****************************************/
    /*              Predictors              */
    /****************************************/
    switch(p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            _primaryPredictor = std::move(bandwidthPredictor);
            _secondaryPredictor = std::move(powerPredictor);
        }break;
        case CONTRACT_POWER_BUDGET:{
            _primaryPredictor = std::move(powerPredictor);
            _secondaryPredictor = std::move(bandwidthPredictor);
        }break;
        default:{
            throw std::runtime_error("Unknown contract.");
        }break;
    }
}

SelectorPredictive::~SelectorPredictive(){
    ;
}

double SelectorPredictive::getPrimaryPrediction(KnobsValues values){
    auto observation = _observedValues.find(getRealValues(_configuration, values));

    // Get real bandwidth from maximum
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            double maxBandwidth = 0.0;
            if(observation != _observedValues.end()){
                maxBandwidth = observation->second.getMaximumBandwidth();
            }else{
                _primaryPredictor->prepareForPredictions();
                maxBandwidth = _primaryPredictor->predict(values, _bandwidthIn->average());
            }
            return _bandwidthIn->average() / maxBandwidth * 100.0;
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            double maxBandwidth = 0.0;
            if(observation != _observedValues.end()){
                maxBandwidth = observation->second.getMaximumBandwidth();
            }else{
                _primaryPredictor->prepareForPredictions();
                maxBandwidth = _primaryPredictor->predict(values, _bandwidthIn->average());
            }
            return min(maxBandwidth, _bandwidthIn->average());
        }break;
        case CONTRACT_POWER_BUDGET:{
            if(observation != _observedValues.end()){
                return observation->second.watts;
            }else{
                _primaryPredictor->prepareForPredictions();
                return _primaryPredictor->predict(values, _bandwidthIn->average());
            }
        }break;
        default:{
            throw std::runtime_error("Unknown contract.");
        }break;
    }
}

double SelectorPredictive::getSecondaryPrediction(KnobsValues values){
    auto observation = _observedValues.find(getRealValues(_configuration, values));

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            if(observation != _observedValues.end()){
                return observation->second.watts;
            }else{
                _secondaryPredictor->prepareForPredictions();
                return _secondaryPredictor->predict(values, _bandwidthIn->average());
            }
        }break;
        case CONTRACT_POWER_BUDGET:{
            double maxBandwidth = 0;
            if(observation != _observedValues.end()){
                maxBandwidth = observation->second.getMaximumBandwidth();
            }else{
                _secondaryPredictor->prepareForPredictions();
                maxBandwidth = _secondaryPredictor->predict(values, _bandwidthIn->average());
            }
            return min(maxBandwidth, _bandwidthIn->average());
        }break;
        default:{
            throw std::runtime_error("Unknown contract.");
        }break;
    }
}

bool SelectorPredictive::isBestSuboptimalValue(double x, double y) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            // Concerning utilization factors, if both are suboptimal,
            // we prefer the closest to the lower bound.
            double distanceX, distanceY;
            distanceX = _p.underloadThresholdFarm - x;
            distanceY = _p.underloadThresholdFarm - y;
            if(distanceX > 0 && distanceY < 0){
                return true;
            }else if(distanceX < 0 && distanceY > 0){
                return false;
            }else{
                return abs(distanceX) < abs(distanceY);
            }
        }break;
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            // Concerning bandwidths, if both are suboptimal,
            // we prefer the higher one.
            return x > y;
        }break;
        case CONTRACT_POWER_BUDGET:{
            // Concerning power budgets, if both are suboptimal,
            // we prefer the lowest one.
            return x < y;
        }break;
        default:{
            ;
        }break;
    }
    return false;
}

bool SelectorPredictive::isBestSecondaryValue(double x, double y) const{
    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_BANDWIDTH:{
            return x < y;
        }break;
        case CONTRACT_POWER_BUDGET:{
            return x > y;
        }break;
        default:{
            ;
        }break;
    }
    return false;
}

KnobsValues SelectorPredictive::getBestKnobsValues(){
    KnobsValues bestValues(KNOB_VALUE_REAL);
    KnobsValues bestSuboptimalValues = _configuration.getRealValues();

    double primaryPrediction = 0;
    double secondaryPrediction = 0;

#ifdef DEBUG_SELECTORS
    double bestPrimaryPrediction = 0;
    double suboptimalSecondary = 0;
#endif
    double bestSecondaryPrediction = 0;
    double bestSuboptimalValue = 0;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            bestSuboptimalValue = _samples->average().bandwidth;
            // We have to minimize the power.
            bestSecondaryPrediction = numeric_limits<double>::max();
        }break;
        case CONTRACT_POWER_BUDGET:{
            bestSuboptimalValue = _samples->average().watts;
            // We have to maximize the bandwidth.
            bestSecondaryPrediction = numeric_limits<double>::min();
        }break;
        default:{
            throw std::runtime_error("Unkown contract type.");
        }break;
    }

    _feasible = false;

    vector<KnobsValues> combinations = _configuration.getAllRealCombinations();
    for(size_t i = 0; i < combinations.size(); i++){
        KnobsValues currentValues = combinations.at(i);
        primaryPrediction = getPrimaryPrediction(currentValues);
        secondaryPrediction = getSecondaryPrediction(currentValues);

        //std::cout << currentValues << " " << primaryPrediction << " ";
        if(isFeasiblePrimaryValue(primaryPrediction, true)){
            //std::cout << secondaryPrediction;
            if(isBestSecondaryValue(secondaryPrediction, bestSecondaryPrediction)){
                bestValues = currentValues;
                _feasible = true;
                DEBUGB(bestPrimaryPrediction = primaryPrediction);
                bestSecondaryPrediction = secondaryPrediction;
            }
        }else if(!_feasible &&
                 isBestSuboptimalValue(primaryPrediction, bestSuboptimalValue)){
            bestSuboptimalValue = primaryPrediction;
            bestSuboptimalValues = currentValues;
            DEBUGB(suboptimalSecondary = secondaryPrediction);
        }
        //std::cout << std::endl;
    }

    if(_feasible){
        DEBUG("Best solution found: " << bestValues);
        DEBUG("Primary prediction: " << bestPrimaryPrediction);
        DEBUG("Secondary prediction: " << bestSecondaryPrediction);
        return bestValues;
    }else{
        DEBUG("Suboptimal solution found: " << bestSuboptimalValues);
        DEBUG("Primary prediction: " << bestSuboptimalValue);
        DEBUG("Secondary prediction: " << suboptimalSecondary);
        return bestSuboptimalValues;
    }
}

void SelectorPredictive::refine(){
    _primaryPredictor->refine();
    _secondaryPredictor->refine();
    _observedValues[_configuration.getRealValues()] = _samples->average();
}

void SelectorPredictive::updatePredictions(const KnobsValues& next){
    KnobsValues real;

    if(next.areRelative()){
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            double realv;
            _configuration.getKnob((KnobType)i)->getRealFromRelative(next[(KnobType)i], realv);
            real[(KnobType)i] = realv;
        }
    }else{
        real = next;
    }

    _primaryPrediction = getPrimaryPrediction(real);
    _secondaryPrediction = getSecondaryPrediction(real);
}

bool SelectorPredictive::predictorsReady() const{
    return _primaryPredictor->readyForPredictions() &&
           _secondaryPredictor->readyForPredictions();
}

bool SelectorPredictive::predictionsDone() const{
    return _primaryPrediction != NOT_VALID &&
           _secondaryPrediction != NOT_VALID;
}

void SelectorPredictive::clearPredictors(){
    _primaryPredictor->clear();
    _secondaryPredictor->clear();
    _primaryPrediction = NOT_VALID;
    _secondaryPrediction = NOT_VALID;
    _observedValues.clear();
}

bool SelectorPredictive::isAccurate(){
    double maxBandwidth = 0.0; double predictedMaxBandwidth = 0.0;
    double power = 0.0; double predictedPower = 0.0;

    maxBandwidth = _samples->average().getMaximumBandwidth();
    power = _samples->average().watts;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            // primaryPrediction is rho, we need to convert to max bandwidth.
            predictedMaxBandwidth = _bandwidthIn->average() / (_primaryPrediction / 100.0);
            predictedPower = _secondaryPrediction;
        }break;
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_BANDWIDTH:{
            predictedMaxBandwidth = _primaryPrediction;
            predictedPower = _secondaryPrediction;
        }break;
        case CONTRACT_POWER_BUDGET:{
            predictedMaxBandwidth = _secondaryPrediction;
            predictedPower = _primaryPrediction;
        }break;
        default:{
            ;
        }break;
    }

    double performanceError = 100.0, powerError = 100.0;
    performanceError = std::abs((maxBandwidth - predictedMaxBandwidth) / maxBandwidth*100.0);
    powerError = std::abs((power - predictedPower) / power*100.0);

    DEBUG("Perf error: " << performanceError);
    DEBUG("Power error: " << powerError);

    if(performanceError > _p.maxPerformancePredictionError ||
       powerError > _p.maxPowerPredictionError /* ||
       _primaryPredictor->getModelError() > 10 || //TODO
       _secondaryPredictor->getModelError() > 10*/){
        return false;
    }else{
        return true;
    }
}

bool SelectorPredictive::isBestSolutionFeasible() const{
    return _feasible;
}

SelectorAnalytical::SelectorAnalytical(const Parameters& p,
               const Configuration& configuration,
               const Smoother<MonitoredSample>* samples):
    SelectorPredictive(p, configuration, samples,
                       std::unique_ptr<Predictor>(new PredictorAnalytical(PREDICTION_BANDWIDTH, p, configuration, samples)),
                       std::unique_ptr<Predictor>(new PredictorAnalytical(PREDICTION_POWER, p, configuration, samples))),
    _violations(0){
    ;
}

KnobsValues SelectorAnalytical::getNextKnobsValues(u_int64_t totalTasks){
    _previousConfiguration = _configuration.getRealValues();
    if(isContractViolated()){
        if(_violations > _p.tolerableSamples){
            _violations = 0;
            return getBestKnobsValues();
        }else{
            ++_violations;
            return _configuration.getRealValues();
        }
    }else{
        if(_violations){
            --_violations;
        }
        return _configuration.getRealValues();
    }
}

std::unique_ptr<Predictor> SelectorLearner::getPredictor(PredictorType type,
                                                         const Parameters& p,
                                                         const Configuration& configuration,
                                                         const Smoother<MonitoredSample>* samples) const{
    Predictor* predictor;
    switch(p.strategyPrediction){
        case STRATEGY_PREDICTION_REGRESSION_LINEAR:{
            predictor = new PredictorLinearRegression(type, p, configuration, samples);
        }break;
        case STRATEGY_PREDICTION_REGRESSION_LINEAR_MAPPING:{
            predictor = new PredictorLinearRegressionMapping(type, p, configuration, samples);
        }break;
        case STRATEGY_PREDICTION_MISHRA:{
            predictor = new PredictorMishra(type, p, configuration, samples);
        }break;
        default:{
            throw std::runtime_error("Unknown prediction strategy.");
        }break;
    }
    return std::unique_ptr<Predictor>(predictor);
}

SelectorLearner::SelectorLearner(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples):
             SelectorPredictive(p, configuration, samples,
                                getPredictor(PREDICTION_BANDWIDTH, p, configuration, samples),
                                getPredictor(PREDICTION_POWER, p, configuration, samples)),
             _explorer(NULL), _firstPointGenerated(false),
             _contractViolations(0), _accuracyViolations(0){
    /***************************************/
    /*              Explorers              */
    /***************************************/
    switch(_p.strategyExploration){
        case STRATEGY_EXPLORATION_RANDOM:{
            _explorer = new ExplorerRandom(configuration);
        }break;
        case STRATEGY_EXPLORATION_HALTON:
        case STRATEGY_EXPLORATION_HALTON_REVERSE:
        case STRATEGY_EXPLORATION_NIEDERREITER:
        case STRATEGY_EXPLORATION_SOBOL:{
            _explorer = new ExplorerLowDiscrepancy(configuration, _p.strategyExploration);
        }break;
        default:{
            throw std::runtime_error("Unknown exploration strategy.");
        }
    }

    if(_p.strategyPrediction == STRATEGY_PREDICTION_REGRESSION_LINEAR_MAPPING){
        _explorer = new ExplorerMultiple(_configuration, _explorer, KNOB_TYPE_MAPPING, MAPPING_TYPE_NUM);
    }
}

SelectorLearner::~SelectorLearner(){
    if(_explorer){
        delete _explorer;
    }
}

KnobsValues SelectorLearner::getNextKnobsValues(u_int64_t totalTasks){
    _previousConfiguration = _configuration.getRealValues();
    KnobsValues kv;
    bool contractViolated = isContractViolated();
    bool accurate = isAccurate();

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
    if(!_firstPointGenerated){
        _firstPointGenerated = true;
        startCalibration(totalTasks);
    }else{
        if(isCalibrating()){
            refine();
        }
    }

    if(isCalibrating()){
        ++_numCalibrationPoints;
        if(!predictorsReady()){
            kv = _explorer->nextRelativeKnobsValues();
        }else{
            if(predictionsDone() && accurate){
                kv = getBestKnobsValues();
                updatePredictions(kv);
                DEBUG("Finished in " << _numCalibrationPoints <<
                      " steps with configuration " << kv);
                stopCalibration(totalTasks);
            }else{
                kv = _explorer->nextRelativeKnobsValues();
                updatePredictions(kv);
            }
        }
    }else{
        if(contractViolated){++_contractViolations;}
        if(!accurate){++_accuracyViolations;}

        if(phaseChanged()){
            /******************* Phase change. *******************/
            _explorer->reset();
            kv = _explorer->nextRelativeKnobsValues();

            // Drop old models.
            clearPredictors();

            startCalibration(totalTasks);
            resetTotalCalibrationTime();
            _accuracyViolations = 0;
            _contractViolations = 0;
            DEBUG("Phase changed, recalibrating.");
        }else if(_bandwidthIn->coefficientVariation() > 100.0){ //TODO Remove magic number
            /******************* Bandwidth change. *******************/
            refine();
            kv = getBestKnobsValues();
            updatePredictions(kv);
            _accuracyViolations = 0;
            _contractViolations = 0;
            _bandwidthIn->reset();
            DEBUG("Input bandwidth fluctuations, recomputing best solution.");
        }else if((!_p.maxCalibrationTime || getTotalCalibrationTime() < _p.maxCalibrationTime) &&
                 ((!accurate && _accuracyViolations > _p.tolerableSamples) ||
                  (contractViolated && _contractViolations > _p.tolerableSamples))){
            /******************* More calibration points. *******************/
            kv = _explorer->nextRelativeKnobsValues();
            updatePredictions(kv);
            refine();
            startCalibration(totalTasks);
            _accuracyViolations = 0;
            _contractViolations = 0;
            if(!accurate){
                DEBUG("Inaccurate model, adding more points.");
            }else{
                DEBUG("Contract violated, adding more points.");
            }
        }else{
            /******************* Stable. *******************/
            if(accurate && _accuracyViolations){ --_accuracyViolations;}
            if(!contractViolated && _contractViolations){ --_contractViolations;}
            kv = _configuration.getRealValues();
        }
    }
    return kv;
}

bool SelectorLearner::phaseChanged() const{
    if(!_configuration.equal(_previousConfiguration)){
        // We need to check that this configuration is equal to the previous one
        // to avoid to detect as a phase change a configuration change.
        return false;
    }
    return _samples->coefficientVariation().latency > 20.0 ||
           _samples->coefficientVariation().watts > 20.0;
}

SelectorFixedExploration::SelectorFixedExploration(const Parameters& p,
               const Configuration& configuration,
               const Smoother<MonitoredSample>* samples,
               std::unique_ptr<Predictor> bandwidthPredictor,
               std::unique_ptr<Predictor> powerPredictor,
               size_t numSamples):
            SelectorPredictive(p, configuration, samples, std::move(bandwidthPredictor), std::move(powerPredictor)){
    const std::vector<KnobsValues>& combinations = _configuration.getAllRealCombinations();
    size_t numConfigurations = combinations.size();
    for(size_t i = 0; i < numConfigurations; i += ceil((double)numConfigurations/numSamples)){
        _confToExplore.push_back(combinations.at(i));
    }
}

SelectorFixedExploration::~SelectorFixedExploration(){
    ;
}

KnobsValues SelectorFixedExploration::getNextKnobsValues(u_int64_t totalTasks){
    _previousConfiguration = _configuration.getRealValues();
    if(_confToExplore.size()){
        if(!isCalibrating()){
            startCalibration(totalTasks);
        }else{
            refine();

        }
        KnobsValues r = _confToExplore.back();
        _confToExplore.pop_back();
        ++_numCalibrationPoints;
        return r;
    }else{
        if(isCalibrating()){
            refine();
            KnobsValues kv = getBestKnobsValues();
            stopCalibration(totalTasks);
            return kv;
        }else{
            return _configuration.getRealValues();
        }
    }
}

SelectorMishra::SelectorMishra(const Parameters& p,
               const Configuration& configuration,
               const Smoother<MonitoredSample>* samples):
        SelectorFixedExploration(p, configuration, samples,
                       std::unique_ptr<Predictor>(new PredictorMishra(PREDICTION_BANDWIDTH, p, configuration, samples)),
                       std::unique_ptr<Predictor>(new PredictorMishra(PREDICTION_POWER, p, configuration, samples)),
                       20){
    ;
}

SelectorMishra::~SelectorMishra(){
    ;
}

SelectorFullSearch::SelectorFullSearch(const Parameters& p,
               const Configuration& configuration,
               const Smoother<MonitoredSample>* samples):
        SelectorFixedExploration(p, configuration, samples,
                       std::unique_ptr<Predictor>(new PredictorFullSearch(PREDICTION_BANDWIDTH, p, configuration, samples)),
                       std::unique_ptr<Predictor>(new PredictorFullSearch(PREDICTION_POWER, p, configuration, samples)),
                       configuration.getAllRealCombinations().size()){
    ;
}

SelectorFullSearch::~SelectorFullSearch(){
    ;
}

Frequency SelectorLiMartinez::findNearestFrequency(Frequency f) const{
    Frequency bestDistance = _availableFrequencies.back();
    Frequency bestFrequency = _availableFrequencies.back();
    for(size_t i = 0; i < _availableFrequencies.size(); i++){
        Frequency distance = std::abs(_availableFrequencies.at(i) - f);
        if(distance < bestDistance){
            bestDistance = distance;
            bestFrequency = _availableFrequencies.at(i);
        }
    }
    return bestFrequency;
}

void SelectorLiMartinez::goRight(){
    _low1 = _low2;
    _high1 = _mid2 - 1;
    _low2 = _mid2 + 1;
    _mid1 = (_low1 + _high1) / 2.0;
    _mid2 = (_low2 + _high2) / 2.0;
}

void SelectorLiMartinez::goLeft(){
    _high1 = _mid1 - 1;
    _low2 = _mid1 + 1;
    _high2 = _high1;
    _mid1 = (_low1 + _high1) / 2.0;
    _mid2 = (_low2 + _high2) / 2.0;
}

SelectorLiMartinez::SelectorLiMartinez(const Parameters& p,
                                           const Configuration& configuration,
                                           const Smoother<MonitoredSample>* samples):
        Selector(p, configuration, samples),
        _firstPointGenerated(false), _low1(0), _mid1(0), _high1(0),
        _low2(0), _mid2(0), _high2(0), _midId(1),
        _availableFrequencies(_p.mammut.getInstanceCpuFreq()->getDomains().back()->getAvailableFrequencies()),
        _currentWatts(DBL_MAX), _optimalWatts(DBL_MAX),
        _optimalFrequency(_availableFrequencies.back()), _optimalWorkers(configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getAllowedValues().back()),
        _currentBw(0),_leftBw(0), _rightBw(0), _improved(false){
    _allowedCores = configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getAllowedValues();
}

SelectorLiMartinez::~SelectorLiMartinez(){
    ;
}

KnobsValues SelectorLiMartinez::getNextKnobsValues(u_int64_t totalTasks){
    _previousConfiguration = _configuration.getRealValues();
    KnobsValues kv(KNOB_VALUE_REAL);

    if(!_firstPointGenerated){
        _firstPointGenerated = true;
        uint maxWorkers = _allowedCores.size(); //_configuration.getKnob(KNOB_TYPE_WORKERS)->getRealValue();
        _low2 = 1;
        _mid2 = maxWorkers / 2.0;
        _high2 = maxWorkers;

        DEBUG("Max workers: " << maxWorkers << " mid2: " << _mid2);
        kv[KNOB_TYPE_VIRTUAL_CORES] = _allowedCores[_mid2 - 1];
        kv[KNOB_TYPE_FREQUENCY] = _availableFrequencies.back();
        _midId = 2;

        startCalibration(totalTasks);

        DEBUG("Generating first point: " << kv);
        ++_numCalibrationPoints;
    }else{
        if(!isCalibrating()){
            return _optimalKv;
        }else if(!isContractViolated()){
            ++_numCalibrationPoints;
            _currentWatts = _samples->average().watts;

            if(_currentWatts < _optimalWatts){
                _improved = true;
                DEBUG("Found a new optimal watts: " << _currentWatts << " vs. " << _optimalWatts);
                _optimalWatts = _currentWatts;
                _optimalFrequency = _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealValue();
                _optimalWorkers = _configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealValue();
                DEBUG("Optimal: " << _optimalWorkers << ", " << _optimalFrequency);
            }

            // We should keep decreasing the frequency
            Frequency currentFrequency = _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealValue();
            kv[KNOB_TYPE_VIRTUAL_CORES] = _configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealValue();

            Frequency nextFrequency = currentFrequency * (_p.requiredBandwidth / _samples->average().bandwidth);
            nextFrequency = findNearestFrequency(nextFrequency);
            if(nextFrequency == currentFrequency){
                --_numCalibrationPoints;
                goto changeworkers;
            }else{
                kv[KNOB_TYPE_FREQUENCY] = nextFrequency;
                DEBUG("Keeping going down on frequencies. We move to: " << kv);
            }
        }else{
changeworkers:
            ++_numCalibrationPoints;
            // I have to change the number of workers
            kv[KNOB_TYPE_FREQUENCY] = _availableFrequencies.back();

            if(_optimalWatts == DBL_MAX){
                // Still I have not found a number of workers that satisfied
                // the time requirement. I increase workers. (Go right).
                kv[KNOB_TYPE_VIRTUAL_CORES] = _allowedCores[_mid2 - 1];
                goRight();
                _midId = 2;
                DEBUG("New interval 1: [" << _low1 << "," << _mid1 << "," << _high1 << "]");
                DEBUG("New interval 2: [" << _low2 << "," << _mid2 << "," << _high2 << "]");
                if(_low1 > _high1 || _low2 > _high2){
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _optimalWorkers;
                    kv[KNOB_TYPE_FREQUENCY] = _optimalFrequency;
                    _optimalKv = kv;
                    stopCalibration(totalTasks);
                    DEBUG("Exploration finished with: " << kv);
                }

            }else if(_currentWatts > _optimalWatts || !_improved){
                DEBUG("This number of workers is worst than the best we found "
                      "up to now.");
                // This number of workers is not ok
                if(_midId == 1){
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _allowedCores[_mid2 - 1];
                    _midId = 2;
                    DEBUG("Trying with the right side. We move to " << kv);
                }else{
                    // Both explored and both are not ok, finished
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _optimalWorkers;
                    kv[KNOB_TYPE_FREQUENCY] = _optimalFrequency;
                    _optimalKv = kv;
                    stopCalibration(totalTasks);
                    DEBUG("Both side are worst. Terminated with: " << kv);
                }
            }else{
                _improved = false;
                if(_midId == 1){
                    goLeft();
                }else{
                    goRight();
                }

                DEBUG("New interval 1: [" << _low1 << "," << _mid1 << "," << _high1 << "]");
                DEBUG("New interval 2: [" << _low2 << "," << _mid2 << "," << _high2 << "]");

                if(_low1 <= _high1){
                    _midId = 1;
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _allowedCores[_mid1 - 1];
                    DEBUG("We move to " << kv);
                }else if(_low2 <= _high2){
                    _midId = 2;
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _allowedCores[_mid2 - 1];
                    DEBUG("We move to " << kv);
                }else{
                    kv[KNOB_TYPE_VIRTUAL_CORES] = _optimalWorkers;
                    kv[KNOB_TYPE_FREQUENCY] = _optimalFrequency;
                    _optimalKv = kv;
                    stopCalibration(totalTasks);
                    DEBUG("Exploration finished with: " << kv);
                }
            }
        }
    }
    return kv;
}

}
