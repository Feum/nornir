/*
 * predictors.cpp
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

/*!
 * \file predictors.cpp
 * \brief Predictors used by adaptive farm.
 */

#include "predictors.hpp"
#include "manager.hpp"

#include "external/Mammut/mammut/cpufreq/cpufreq.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_PREDICTORS
#define DEBUG(x) do { std::cerr << "[Predictors] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif


namespace nornir{

using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

using namespace mlpack::regression;

using namespace arma;

static double getVoltage(VoltageTable table, uint workers, Frequency frequency){
    VoltageTableKey key(workers, frequency);
    VoltageTableIterator it = table.find(key);
    if(it != table.end()){
        return it->second;
    }else{
        throw runtime_error("Frequency and/or number of virtual cores "
                            "not found in voltage table.");
    }
}

/**
 * Computes a proportional value for power consumption of a configuration.
 * ATTENTION: Assumes linear mapping of the workers on the cores.
 */
static void getPowerProportions(VoltageTable table, uint physicalCores,
                                Frequency frequency, uint coresPerDomain,
                                double& staticPower, double& dynamicPower){
    int numCores = physicalCores;
    double voltage = 0;
    staticPower = 0;
    dynamicPower = 0;
    uint currentCores = 0;
    while(numCores > 0){
        currentCores = (numCores < (int) coresPerDomain)?(uint)numCores:coresPerDomain;
        voltage = getVoltage(table, currentCores, frequency);

        staticPower += voltage;
        dynamicPower += currentCores*frequency*voltage*voltage;

        numCores -= coresPerDomain;
    }
}

RegressionData::RegressionData(const Parameters& p,
                               const FarmConfiguration& configuration,
                               const Smoother<MonitoredSample>* samples):
        _p(p), _configuration(configuration), _samples(samples){
    _topology = _p.mammut.getInstanceTopology();
    _cpus = _topology->getCpus().size();
    _phyCores = _topology->getPhysicalCores().size();
    _phyCoresPerCpu = _topology->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = _topology->getPhysicalCore(0)->getVirtualCores().size();
}

double RegressionData::getUsedPhysicalCores(double numWorkers, bool includeServiceNodes){
    uint numNodes = numWorkers;
    if(includeServiceNodes){
        numNodes += _configuration.getNumServiceNodes();
    }
    return std::min(numNodes, _phyCores);
}

void RegressionData::init(){init(_configuration.getRealValues());}

void RegressionDataServiceTime::init(const KnobsValues& values){
    _numPredictors = 0;
    uint workersCores = values[KNOB_TYPE_WORKERS];
    double physicalCores = getUsedPhysicalCores(workersCores, false);

    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        double frequency = values[KNOB_TYPE_FREQUENCY];
        _invScalFactorFreq = (double)_minFrequency / frequency;
        ++_numPredictors;

        _invScalFactorFreqAndCores = (double)_minFrequency /
                                     (physicalCores * frequency);
        ++_numPredictors;
    }
}

RegressionDataServiceTime::RegressionDataServiceTime(const Parameters& p,
                                                     const FarmConfiguration& configuration,
                                                     const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples),
        _invScalFactorFreq(0),
        _invScalFactorFreqAndCores(0),
        _numPredictors(0){
    _phyCores = _p.mammut.getInstanceTopology()->getPhysicalCores().size();
    _minFrequency = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->
                                                    getAvailableFrequencies().at(0);
    init(_configuration.getRealValues());
}

uint RegressionDataServiceTime::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataServiceTime::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        matrix(rowId++, columnId) = _invScalFactorFreq;
        matrix(rowId++, columnId) = _invScalFactorFreqAndCores;
    }
}

void RegressionDataPower::init(const KnobsValues& values){
    _numPredictors = 0;
    Frequency frequency = values[KNOB_TYPE_FREQUENCY];
    uint workersCores = values[KNOB_TYPE_WORKERS];
    double usedPhysicalCores = getUsedPhysicalCores(workersCores, true);

    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        uint usedCpus = 0;
        if(_p.knobMapping == KNOB_MAPPING_LINEAR){
            switch(_p.knobHyperthreading){
                case KNOB_HT_AUTO:{
                    throw std::runtime_error("This should never happen.");
                }
                case KNOB_HT_LATER:
                case KNOB_HT_SOONER:
                case KNOB_HT_NO:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         (double) _phyCoresPerCpu);
                }break;
                /*
                case KNOB_HT_LATER:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         (double) _phyCoresPerCpu);
                    if(usedPhysicalCores > _phyCores){
                        _additionalContextes = usedPhysicalCores - _phyCores;
                        _additionalContextes = std::min(_additionalContextes, _virtCoresPerPhyCores*_phyCores);
                    }else{
                        _additionalContextes = 0;
                    }
                    ++_numPredictors;
                }break;
                case KNOB_HT_SOONER:{
                    usedCpus = std::ceil(usedPhysicalCores /
                                         ((double) _phyCoresPerCpu * (double) _virtCoresPerPhyCores));
                }break;
                */
            }
        }else{
            ; //TODO
        }
        usedCpus = std::min(usedCpus, _cpus);
        uint unusedCpus = _cpus - usedCpus;

        double staticPowerProp = 0, dynamicPowerProp = 0;
        //TODO: Potrebbe non esserci una frequenza se FREQUENCY_NO. In tal caso non possiamo predirre nulla.
        getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                frequency, _phyCoresPerCpu,
                staticPowerProp, dynamicPowerProp);
        _voltagePerUsedSockets = staticPowerProp;
        ++_numPredictors;

        if(_cpus > 1){
            Frequency frequencyUnused;
            switch(_strategyUnused){
                /**
                 * TODO:
                 * For simplicity we assume that the frequency of the unused
                 * sockets is equal to the frequency of the used sockets. This
                 * should be modified as future work.
                 */
                default:{
                    frequencyUnused = frequency;
                }
            }
            double voltage = getVoltage(_p.archData.voltageTable, 0, frequencyUnused);
            _voltagePerUnusedSockets = voltage * unusedCpus;
            ++_numPredictors;
        }

        _dynamicPowerModel = dynamicPowerProp;
        ++_numPredictors;
    }else{
        /**
         * Since I do not control the frequency, this model can be very
         * inaccurate.
         */
        _dynamicPowerModel = usedPhysicalCores;
        ++_numPredictors;
    }
}

RegressionDataPower::RegressionDataPower(const Parameters& p,
                                         const FarmConfiguration& configuration,
                                         const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples), _dynamicPowerModel(0),
        _voltagePerUsedSockets(0), _voltagePerUnusedSockets(0),
        _additionalContextes(0), _numPredictors(0){
    _topology = _p.mammut.getInstanceTopology();
    _cpus = _topology->getCpus().size();
    _phyCores = _topology->getPhysicalCores().size();
    _phyCoresPerCpu = _topology->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = _topology->getPhysicalCore(0)->getVirtualCores().size();
    _strategyUnused = _p.strategyUnusedVirtualCores;
    init(_configuration.getRealValues());
}

uint RegressionDataPower::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataPower::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    matrix(rowId++, columnId) = _dynamicPowerModel;
    if(_p.knobFrequencies == KNOB_FREQUENCY_YES){
        matrix(rowId++, columnId) = _voltagePerUsedSockets;
        if(_cpus > 1){
            matrix(rowId++, columnId) = _voltagePerUnusedSockets;
        }
    }
    /*
    if(_p.knobHyperthreading != KNOB_HT_NO){
        matrix(rowId++, columnId) = _additionalContextes;
    }
    */
}

Predictor::Predictor(PredictorType type,
                     const Parameters& p,
                     const FarmConfiguration& configuration,
                     const Smoother<MonitoredSample>* samples):
        _type(type), _p(p), _configuration(configuration),
        _samples(samples), _modelError(0){
    ;
}

Predictor::~Predictor(){
    ;
}

PredictorLinearRegression::PredictorLinearRegression(PredictorType type,
                                                     const Parameters& p,
                                                     const FarmConfiguration& configuration,
                                                     const Smoother<MonitoredSample>* samples):
        Predictor(type, p, configuration, samples), _preparationNeeded(true){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            _predictionInput = new RegressionDataServiceTime(p, configuration, samples);
        }break;
        case PREDICTION_POWER:{
            _predictionInput = new RegressionDataPower(p, configuration, samples);
        }break;
    }
    clear();
}

PredictorLinearRegression::~PredictorLinearRegression(){
    clear();
    delete _predictionInput;
}

void PredictorLinearRegression::clear(){
    for(obs_it iterator = _observations.begin();
               iterator != _observations.end();
               iterator++){
        delete iterator->second.data;
    }
    _observations.clear();
    _agingVector.clear();
    _agingVector.reserve(_p.regressionAging);
    _currentAgingId = 0;
}

bool PredictorLinearRegression::readyForPredictions(){
    _predictionInput->init(_configuration.getRealValues());
    uint minPoints = std::max<uint>(_predictionInput->getNumPredictors(), 2);
    return _observations.size() >= minPoints;
}

double PredictorLinearRegression::getCurrentResponse() const{
    double r = 0.0;
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            r = 1.0 / _samples->average().bandwidth;
        }break;
        case PREDICTION_POWER:{
            r = _samples->average().watts;
        }break;
    }
    return r;
}

void PredictorLinearRegression::refine(){
    KnobsValues currentValues = _configuration.getRealValues();
    _preparationNeeded = true;
    if(_type == PREDICTION_POWER){
        // Add service nodes
        currentValues[KNOB_TYPE_WORKERS] = currentValues[KNOB_TYPE_WORKERS] + _configuration.getNumServiceNodes();
    }
    obs_it lb = _observations.lower_bound(currentValues);

    if(_p.regressionAging && !contains(_agingVector, currentValues)){
        if(_agingVector.size() < _p.regressionAging){
            _agingVector.push_back(currentValues);
        }else{
            _agingVector.at(_currentAgingId) = currentValues;
        }
        _currentAgingId = (_currentAgingId + 1) % _p.regressionAging;
    }
    DEBUG("Refining with configuration " << currentValues << ": "
                                         << getCurrentResponse());
    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentValues, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentValues);
        lb->second.data->init();
        lb->second.response = getCurrentResponse();
    }else{
        // The key does not exist in the map
        Observation o;
        switch(_type){
            case PREDICTION_BANDWIDTH:{
                o.data = new RegressionDataServiceTime(_p, _configuration, _samples);
            }break;
            case PREDICTION_POWER:{
                o.data = new RegressionDataPower(_p, _configuration, _samples);
            }break;
        }
        o.response = getCurrentResponse();
        _observations.insert(lb, Observations::value_type(currentValues, o));
    }
}

void PredictorLinearRegression::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("prepareForPredictions: Not enough "
                                     "points are present");
        }

        // One observation per column.
        arma::mat dataMl(_observations.begin()->second.data->getNumPredictors(),
                         _observations.size());
        arma::vec responsesMl(_observations.size());

        size_t i = 0;
        for(obs_it iterator = _observations.begin();
                   iterator != _observations.end();
                   iterator++){
            const Observation& obs = iterator->second;

            if(!_p.regressionAging || contains(_agingVector, iterator->first)){
                obs.data->toArmaRow(i, dataMl);
                responsesMl(i) = obs.response;
                ++i;
            }
        }

        if(_p.regressionAging && _agingVector.size() != _observations.size()){
            dataMl.resize(_observations.begin()->second.data->getNumPredictors(),
                          _agingVector.size());
            responsesMl.resize(_agingVector.size());
        }

        _lr = LinearRegression(dataMl, responsesMl);
        _modelError =  _lr.ComputeError(dataMl, responsesMl);
        DEBUG("Error in model: " << _modelError);
        _preparationNeeded = false;
    }
}

double PredictorLinearRegression::predict(const KnobsValues& values){
    _predictionInput->init(values);

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

PredictorAnalytical::PredictorAnalytical(PredictorType type,
                                 const Parameters& p,
                                 const FarmConfiguration& configuration,
                                 const Smoother<MonitoredSample>* samples):
            Predictor(type, p, configuration, samples) {
    Topology* t = _p.mammut.getInstanceTopology();
    _phyCores = t->getPhysicalCores().size();
    _phyCoresPerCpu = t->getCpu(0)->getPhysicalCores().size();
}

double PredictorAnalytical::getScalingFactor(const KnobsValues& values){
    double usedPhysicalCores = values[KNOB_TYPE_WORKERS];
    return (double)(values[KNOB_TYPE_FREQUENCY] * usedPhysicalCores) /
           (double)(_configuration.getRealValue(KNOB_TYPE_FREQUENCY) *
                    _configuration.getRealValue(KNOB_TYPE_WORKERS));
}

double PredictorAnalytical::getPowerPrediction(const KnobsValues& values){
    double staticPower = 0, dynamicPower = 0;
    double usedPhysicalCores = values[KNOB_TYPE_WORKERS];
    getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                        values[KNOB_TYPE_FREQUENCY], _phyCoresPerCpu,
                        staticPower, dynamicPower);
    return dynamicPower;
}

void PredictorAnalytical::prepareForPredictions(){
    ;
}

bool PredictorAnalytical::readyForPredictions(){
    return true;
}

void PredictorAnalytical::refine(){
    ;
}

void PredictorAnalytical::clear(){
    ;
}

double PredictorAnalytical::predict(const KnobsValues& values){
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            return _samples->average().bandwidth * getScalingFactor(values);
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(values);
        }break;
    }
    return 0.0;
}


PredictorMishra::PredictorMishra(PredictorType type,
              const Parameters& p,
              const FarmConfiguration& configuration,
              const Smoother<MonitoredSample>* samples):
                  Predictor(type, p, configuration, samples),
                  _preparationNeeded(true){
    const std::vector<KnobsValues>& combinations = _configuration.getAllRealCombinations();
    DEBUG("Found: " << combinations.size() << " combinations.");
    for(size_t i = 0; i < combinations.size(); i++){
        _confIndexes[combinations.at(i)] = i;
    }
    _values.resize(combinations.size());
    _values.zeros();
}

PredictorMishra::~PredictorMishra(){
    ;
}

bool PredictorMishra::readyForPredictions(){
    return true;
}

void PredictorMishra::clear(){
    _values.zeros();
    _predictions.zeros();
}

void PredictorMishra::refine(){
    _preparationNeeded = true;
    auto it = _confIndexes.find(_configuration.getRealValues());
    if(it == _confIndexes.end()){
        throw std::runtime_error("[Mishra] Impossible to find index for configuration.");
    }
    size_t confId = it->second;
    if(confId >= _values.size()){
        throw std::runtime_error("[Mishra] Invalid configuration index: " + confId);
    }
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            _values.at(confId) = _samples->average().bandwidth;
        }break;
        case PREDICTION_POWER:{
            _values.at(confId) = _samples->average().watts;
        }break;
        default:{
            throw std::runtime_error("[Mishra] Unknown predictor type.");
        }
    }
}

void PredictorMishra::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("[Mishra] prepareForPredictions: Not enough "
                                     "points are present");
        }

        string dataFile;
        bool perColumnNormalization;
        switch(_type){
            case PREDICTION_BANDWIDTH:{
                dataFile = _p.mishra.bandwidthData;
                perColumnNormalization = true;
            }break;
            case PREDICTION_POWER:{
                dataFile = _p.mishra.powerData;
                perColumnNormalization = false;
            }break;
            default:{
                throw std::runtime_error("Unknown predictor type.");
            }
        }

        leo::PredictionResults pr = leo::compute(_p.mishra.appId, dataFile,
                                                 &_values, perColumnNormalization);
        _predictions = arma::conv_to<std::vector<double> >::from(pr.predictions);
#if 0
        for(arma::vec::iterator i = pr.predictions.begin(); i != pr.predictions.end(); ++i){
            _predictions.push_back(*i);
        }
#endif
        _preparationNeeded = false;
    }
}

double PredictorMishra::predict(const KnobsValues& realValues){
    auto it = _confIndexes.find(realValues);
    if(it == _confIndexes.end()){
        throw std::runtime_error("[Mishra] Impossible to find index for configuration.");
    }
    size_t confId = it->second;
    if(confId >= _predictions.size()){
        throw std::runtime_error("[Mishra] Invalid configuration index: " + confId);
    }
    return _predictions.at(confId);
}

PredictorFullSearch::PredictorFullSearch(PredictorType type,
          const Parameters& p,
          const FarmConfiguration& configuration,
          const Smoother<MonitoredSample>* samples):
      Predictor(type, p, configuration, samples),
      _allConfigurations(_configuration.getAllRealCombinations()){
    ;
}

PredictorFullSearch::~PredictorFullSearch(){
    ;
}

bool PredictorFullSearch::readyForPredictions(){
    return _values.size() == _allConfigurations.size();
}

void PredictorFullSearch::clear(){
    _values.clear();
}

void PredictorFullSearch::refine(){
    double value = 0;
    switch(_type){
        case PREDICTION_BANDWIDTH:{
            value = _samples->average().bandwidth;
        }break;
        case PREDICTION_POWER:{
            value = _samples->average().watts;
        }break;
        default:{
            throw std::runtime_error("Unknown predictor type.");
        }
    }
    _values[_configuration.getRealValues()] = value;
}

void PredictorFullSearch::prepareForPredictions(){
    ;
}

double PredictorFullSearch::predict(const KnobsValues& realValues){
    if(!readyForPredictions()){
        throw std::runtime_error("prepareForPredictions: Not enough "
                                 "points are present");
    }
    return _values.at(realValues);
}

}

