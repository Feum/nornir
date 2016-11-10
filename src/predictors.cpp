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
 */
static void getPowerProportions(VoltageTable table, uint physicalCores,
                                Frequency frequency, uint coresPerDomain,
                                uint numDomains, MappingType mt,
                                double& staticPower, double& dynamicPower){
    int numCores = physicalCores;
    double voltage = 0;
    staticPower = 0;
    dynamicPower = 0;
    uint currentCores = 0;
    int interleavedPerDomain = physicalCores / numDomains;
    int interleavedSpurious = physicalCores % numDomains;
    while(numCores > 0){
        switch(mt){
            case MAPPING_TYPE_LINEAR:{
                currentCores = (numCores < (int) coresPerDomain)?(uint)numCores:coresPerDomain;
            }break;
            case MAPPING_TYPE_INTERLEAVED:{
                currentCores =  interleavedPerDomain;
                if(interleavedSpurious){
                    currentCores +=1;
                    interleavedSpurious--;
                }
            }break;
            default:{
                throw std::runtime_error("getPowerProportions: mapping type not supported.");
            }
        }
        voltage = getVoltage(table, currentCores, frequency);

        staticPower += voltage;
        dynamicPower += currentCores*frequency*voltage*voltage;

        numCores -= currentCores;
    }
}

RegressionData::RegressionData(const Parameters& p,
                               const Configuration& configuration,
                               const Smoother<MonitoredSample>* samples):
        _p(p), _configuration(configuration), _samples(samples){
    _topology = _p.mammut.getInstanceTopology();
    _cpus = _topology->getCpus().size();
    _phyCores = _topology->getPhysicalCores().size();
    _phyCoresPerCpu = _topology->getCpu(0)->getPhysicalCores().size();
    _virtCoresPerPhyCores = _topology->getPhysicalCore(0)->getVirtualCores().size();
}

double RegressionData::getUsedPhysicalCores(double numVirtualCores){
    if(_configuration.getRealValue(KNOB_TYPE_HYPERTHREADING) != 1){
        throw std::runtime_error("getUsedPhysicalCores, ht > 1 not yet implemented.");
    }
    uint realVirtualCores = numVirtualCores + _configuration.getNumServiceNodes();
    return std::min(realVirtualCores, _phyCores);
}

void RegressionData::init(){init(_configuration.getRealValues());}

void RegressionDataServiceTime::init(const KnobsValues& values){
    _numPredictors = 0;
    uint numVirtualCores = values[KNOB_TYPE_VIRTUAL_CORES];

    if(_p.knobFrequencyEnabled){
        double frequency = values[KNOB_TYPE_FREQUENCY];
        _invScalFactorFreq = (double)_minFrequency / frequency;
        ++_numPredictors;

        _invScalFactorFreqAndCores = (double)_minFrequency /
                                     (numVirtualCores * frequency);
        ++_numPredictors;
    }else{
        throw std::runtime_error("Impossible to use Amdahl regression with no frequency.");
    }
}

RegressionDataServiceTime::RegressionDataServiceTime(const Parameters& p,
                                                     const Configuration& configuration,
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
    if(_p.knobFrequencyEnabled){
        matrix(rowId++, columnId) = _invScalFactorFreq;
        matrix(rowId++, columnId) = _invScalFactorFreqAndCores;
    }
}

void RegressionDataPower::init(const KnobsValues& values){
    assert(values.areReal());
    _numPredictors = 0;
    Frequency frequency = values[KNOB_TYPE_FREQUENCY];
    uint numVirtualCores = values[KNOB_TYPE_VIRTUAL_CORES];
    double usedPhysicalCores = getUsedPhysicalCores(numVirtualCores);

    if(_p.knobFrequencyEnabled){
        uint usedCpus = 0;
        if(values[KNOB_TYPE_MAPPING] == MAPPING_TYPE_LINEAR){
            usedCpus = std::ceil(usedPhysicalCores / (double) _phyCoresPerCpu);
        }else{
            usedCpus = std::min(usedPhysicalCores, (double) _cpus);
        }
        usedCpus = std::min(usedCpus, _cpus);
        uint unusedCpus = _cpus - usedCpus;

        double staticPowerProp = 0, dynamicPowerProp = 0;
        //TODO: Potrebbe non esserci una frequenza se FREQUENCY_NO. In tal caso non possiamo predirre nulla.
        getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                frequency, _phyCoresPerCpu,
                _cpus, (MappingType) values[KNOB_TYPE_MAPPING],
                staticPowerProp, dynamicPowerProp);
        //_voltagePerUsedSockets = staticPowerProp;
        //++_numPredictors;

        if(_cpus > 1){
            Frequency frequencyUnused;
            switch(_strategyUnused){
                case STRATEGY_UNUSED_VC_SAME:{
                    frequencyUnused = frequency;
                }break;
                case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                    frequencyUnused = _lowestFrequency;
                }break;
                default:{
                    throw std::runtime_error("RegressionDataPower: init: strategyUnused unsupported.");
                }
            }
            double voltage = getVoltage(_p.archData.voltageTable, 0, frequencyUnused);
            _voltagePerUnusedSockets = (voltage * unusedCpus) + staticPowerProp;
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
                                         const Configuration& configuration,
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
    _lowestFrequency = _p.mammut.getInstanceCpuFreq()->getDomains().front()->getAvailableFrequencies().front();
    init(_configuration.getRealValues());
}

uint RegressionDataPower::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataPower::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    matrix(rowId++, columnId) = _dynamicPowerModel;
    if(_p.knobFrequencyEnabled){
        //matrix(rowId++, columnId) = _voltagePerUsedSockets;
        if(_cpus > 1){
            matrix(rowId++, columnId) = _voltagePerUnusedSockets;
        }
    }
}

bool RegressionDataPower::getInactivePowerPosition(size_t& pos) const{
    // ATTENTION: pos = 0 is the intercept (constant value).
    if(_p.knobFrequencyEnabled && _cpus > 1){
        pos = _numPredictors;
        return true;
    }else{
        return false;
    }
}

Predictor::Predictor(PredictorType type,
                     const Parameters& p,
                     const Configuration& configuration,
                     const Smoother<MonitoredSample>* samples):
        _type(type), _p(p), _configuration(configuration),
        _samples(samples), _modelError(0){
    ;
}

Predictor::~Predictor(){
    ;
}

double Predictor::getMaximumBandwidth() const{
    if(_samples->average().utilisation >= MAX_RHO){
        return _samples->average().bandwidth;
    }else{
        return _samples->average().getMaximumBandwidth();
    }
}

double Predictor::getCurrentPower() const{
    return _samples->average().watts;
}

PredictorLinearRegression::PredictorLinearRegression(PredictorType type,
                                                     const Parameters& p,
                                                     const Configuration& configuration,
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
            r = 1.0 / getMaximumBandwidth();
        }break;
        case PREDICTION_POWER:{
            r = getCurrentPower();
        }break;
    }
    return r;
}

void PredictorLinearRegression::refine(){
    KnobsValues currentValues = _configuration.getRealValues();
    _preparationNeeded = true;
    if(_type == PREDICTION_POWER){
        // Add service nodes
        currentValues[KNOB_TYPE_VIRTUAL_CORES] += _configuration.getNumServiceNodes();
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
    double response = getCurrentResponse();
    DEBUG("Refining with configuration " << currentValues << ": "
                                         << response);
    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentValues, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentValues);
        lb->second.data->init();
        lb->second.response = response;
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
        o.response = response;
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

        // Remove rows in dataMl if they contain all the same value.
        // In this way we ensure that the matrix can be inverted and that
        // the regression has a solution.
        size_t id = 0;
        _removedRows.clear();
        for(size_t i = 0; i < dataMl.n_rows; ){
            double firstElem = dataMl(i, 0);
            bool removeRow = true;
            for(size_t j = 1; j < dataMl.n_cols; j++){
                if(dataMl(i, j) != firstElem){
                    removeRow = false;
                    break;
                }
            }
            if(removeRow){
                dataMl.shed_row(i);
                _removedRows.push_back(id);
            }else{
                i++;
            }
            ++id;
        }
        //// End of row removal ////

        if(_p.regressionAging && _agingVector.size() != _observations.size()){
            dataMl.resize(_observations.begin()->second.data->getNumPredictors(),
                          _agingVector.size());
            responsesMl.resize(_agingVector.size());
        }

        std::ostringstream x1;
        std::ostringstream x2;
        arma::set_stream_err1(x1);
        arma::set_stream_err2(x2);
        LinearRegression newlr = LinearRegression(dataMl, responsesMl);
        if(x1.str().compare("") || x2.str().compare("")){
            ;
        }else{
            _lr = newlr;
            _modelError =  _lr.ComputeError(dataMl, responsesMl);
        }

        DEBUG("Error in model: " << _modelError);
        _preparationNeeded = false;
    }
}

double PredictorLinearRegression::predict(const KnobsValues& values){
    _predictionInput->init(values);

    // One observation per column.
    arma::mat predictionInputMl(_predictionInput->getNumPredictors(), 1);
    arma::vec result(1);
    double res = 0;
    _predictionInput->toArmaRow(0, predictionInputMl);

    // Remove from prediction input the corresponding positions
    //for(size_t i = 0; i < _removedRows.size(); i++){
    //    predictionInputMl.shed_row(_removedRows.at(i));
    //}
    arma::mat predictionInputMlShed(_predictionInput->getNumPredictors(), 1);
    size_t realSize = 0;
    for(size_t i = 0; i < _predictionInput->getNumPredictors(); i++){
        if(!contains(_removedRows, i)){
            predictionInputMlShed(realSize, 0) = predictionInputMl(i, 0);
            ++realSize;
        }
    }
    predictionInputMlShed.resize(realSize, 1);

    _lr.Predict(predictionInputMlShed, result);
    if(_type == PREDICTION_BANDWIDTH){
        res = 1.0 / result.at(0);
    }else{
        res = result.at(0);
    }

    return res;
}

double PredictorLinearRegression::getInactivePowerParameter() const{
    if(_type == PREDICTION_POWER){
        size_t pos;
        if(!((RegressionDataPower*)_predictionInput)->getInactivePowerPosition(pos)){
            throw std::runtime_error("Impossible to get inactive power parameter.");
        }
        return _lr.Parameters().at(pos);
    }else{
        throw std::runtime_error("getInactivePowerParameter can only be called on POWER predictors.");
    }
}

/**************** PredictorUSL ****************/
PredictorUsl::PredictorUsl(PredictorType type,
             const Parameters& p,
             const Configuration& configuration,
             const Smoother<MonitoredSample>* samples):
                    Predictor(type, p, configuration, samples),
                    _ws(NULL), _x(NULL), _y(NULL), _chisq(0),
                    _preparationNeeded(true), _maxFreqBw(0), _minFreqBw(0),
					_minFreqCoresBw(0){
    if(type != PREDICTION_BANDWIDTH){
        throw std::runtime_error("PredictorUsl can only be used for bandwidth predictions.");
    }
    _c = gsl_vector_alloc(POLYNOMIAL_DEGREE_USL);
    _cov = gsl_matrix_alloc(POLYNOMIAL_DEGREE_USL, POLYNOMIAL_DEGREE_USL);
    _minFrequency = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->getAvailableFrequencies().front();
    _maxFrequency = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->getAvailableFrequencies().back();
    _maxCores = _configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getAllowedValues().back();
}

PredictorUsl::~PredictorUsl(){
    gsl_matrix_free(_cov);
    gsl_vector_free(_c);
}

void PredictorUsl::clear(){
    _xs.clear();
    _ys.clear();
}

bool PredictorUsl::readyForPredictions(){
    return _xs.size() >= POLYNOMIAL_DEGREE_USL;
}

void PredictorUsl::refine(){
    double numCores = _configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealValue();
    double frequency = _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealValue();
    double bandwidth = getMaximumBandwidth();


    if(frequency == _maxFrequency && numCores == _maxCores){
        _maxFreqBw = bandwidth;
        return;
    }else if(frequency == _minFrequency){
    	if(numCores == _maxCores){
    		_minFreqBw = bandwidth;
    	}else if(numCores == 1){
    		_minFreqCoresBw = bandwidth;
    	}
    }else if(frequency != _minFrequency){
        return;
    }

    double x, y;
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
    	x = numCores - 1;
    	y = ((_minFreqCoresBw * numCores)/bandwidth) - 1;
    }else{
    	x = numCores - 1;
    	y = numCores / bandwidth;
    }

    // Checks if a y is already present for this x
    int pos = -1;
    for(size_t i = 0; i < _xs.size(); i++){
        if(_xs[i] == x){
            pos = i;
            break;
        }
    }
    if(pos == -1){
        _xs.push_back(x);
        _ys.push_back(y);
    }else{
        _ys.at(pos) = y;
    }
    _preparationNeeded = true;
}

void PredictorUsl::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("PredictorUsl: Not yet ready for predictions.");
        }

        uint maxPolDegree = POLYNOMIAL_DEGREE_USL;
        if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        	maxPolDegree -= 1; //To remove x^0 value.
        }

        _x = gsl_matrix_alloc(_xs.size(), maxPolDegree);
        _y = gsl_vector_alloc(_xs.size());
        size_t i = 0, j = 0;
        for(i = 0; i < _xs.size(); i++) {
            for(j = 0; j < POLYNOMIAL_DEGREE_USL; j++) {
            	if(j == 0 && _p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
            		continue; //To skip x^0 value.
            	}
                gsl_matrix_set(_x, i, j, pow(_xs.at(i), j));
            }
            gsl_vector_set(_y, i, _ys.at(i));
        }

        _ws = gsl_multifit_linear_alloc(_xs.size(), maxPolDegree);
        gsl_multifit_linear(_x, _y, _c, _cov, &_chisq, _ws);

        _coefficients.clear();
        for(i = 0; i < maxPolDegree; i++){
            _coefficients.push_back(gsl_vector_get(_c, i));
        }
        gsl_matrix_free(_x);
        gsl_vector_free(_y);
        gsl_multifit_linear_free(_ws);
        _preparationNeeded = false;
    }else{
    	if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
    		std::cout << "B1Pred(actual): " << _minFreqCoresBw << std::endl;
            std::cout << "Contention: " << _coefficients[0] - _coefficients[1] << std::endl;
            std::cout << "Coherency: " << (double) _coefficients[1]/(_coefficients[0] - _coefficients[1]) << std::endl;
    	}else{
            std::cout << "B1Pred: " << 1.0 / _coefficients[0] << std::endl;
            std::cout << "Contention: " << (_coefficients[1] - _coefficients[2])/((double) _coefficients[0]) << std::endl;
            std::cout << "Coherency: " << (double) _coefficients[2]/(_coefficients[1] - _coefficients[2]) << std::endl;
    	}
    }
}

double PredictorUsl::predict(const KnobsValues& configuration){
    double result = 0;
    double numCores = 0;
    double frequency = 0;
    if(configuration.areReal()){
        numCores = configuration[KNOB_TYPE_VIRTUAL_CORES];
        frequency = configuration[KNOB_TYPE_FREQUENCY];
    }else{
        _configuration.getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealFromRelative(configuration[KNOB_TYPE_VIRTUAL_CORES], numCores);
        _configuration.getKnob(KNOB_TYPE_FREQUENCY)->getRealFromRelative(configuration[KNOB_TYPE_FREQUENCY], frequency);
    }
    for(size_t i = 0; i < _coefficients.size(); i++){
        result += _coefficients.at(i)*std::pow(numCores - 1, i);
    }

    double bandwidth;
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
    	bandwidth = (numCores * _minFreqCoresBw)/(result + 1);
    }else{
    	bandwidth = (numCores / result);
    }

    if(frequency != _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->getAvailableFrequencies().front()){
        double minMaxScaling = _maxFreqBw / _minFreqBw;
        double maxFreqPred = bandwidth * minMaxScaling;
        return ((maxFreqPred - bandwidth)/(_maxFrequency - _minFrequency))*frequency + ((bandwidth*_maxFrequency)-(maxFreqPred*_minFrequency))/(_maxFrequency - _minFrequency);
    }else{
        return bandwidth;
    }
}
    
/**************** PredictorSimple ****************/

PredictorAnalytical::PredictorAnalytical(PredictorType type,
                                 const Parameters& p,
                                 const Configuration& configuration,
                                 const Smoother<MonitoredSample>* samples):
            Predictor(type, p, configuration, samples) {
    Topology* t = _p.mammut.getInstanceTopology();
    _phyCores = t->getPhysicalCores().size();
    _phyCoresPerCpu = t->getCpu(0)->getPhysicalCores().size();
    _cpus = t->getCpus().size();
}

double PredictorAnalytical::getScalingFactor(const KnobsValues& values){
    double usedVirtualCores = values[KNOB_TYPE_VIRTUAL_CORES];
    return (double)(values[KNOB_TYPE_FREQUENCY] * usedVirtualCores) /
           (double)(_configuration.getRealValue(KNOB_TYPE_FREQUENCY) *
                    _configuration.getRealValue(KNOB_TYPE_VIRTUAL_CORES));
}

double PredictorAnalytical::getPowerPrediction(const KnobsValues& values){
    assert(values.areReal());
    double staticPower = 0, dynamicPower = 0;
    double usedPhysicalCores = values[KNOB_TYPE_VIRTUAL_CORES];
    getPowerProportions(_p.archData.voltageTable, usedPhysicalCores,
                        values[KNOB_TYPE_FREQUENCY], _phyCoresPerCpu,
                        _cpus, (MappingType) values[KNOB_TYPE_MAPPING],
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
            double scalingFactor = getScalingFactor(values);
            return getMaximumBandwidth() * scalingFactor;
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(values);
        }break;
    }
    return 0.0;
}


PredictorMishra::PredictorMishra(PredictorType type,
              const Parameters& p,
              const Configuration& configuration,
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
    std::vector<std::string> names = mammut::utils::readFile(p.mishra.namesData);
    bool appFound = false;
    for(size_t i = 0; i < names.size(); i++){
        if(names.at(i).compare(p.mishra.applicationName) == 0){
            _appId = i;
            appFound = true;
        }
    }
    if(!appFound){
        throw std::runtime_error("Impossible to find application " + p.mishra.applicationName +
                                 " in " + p.mishra.namesData);
    }
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
            // TODO In the offline profiling data we have bandwidthMax, not bandwidth
            // According, if we do not have PERF_UTILIZATION contracts and we have
            // a utilization < 1, results may be wrong.
#if 0
            if(_samples->average().utilisation < MAX_RHO){
                throw std::runtime_error("[Mishra] bandwidth bandwidthMax problem.");
            }
#endif
            _values.at(confId) = getMaximumBandwidth();
        }break;
        case PREDICTION_POWER:{
            _values.at(confId) = getCurrentPower();
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

        leo::PredictionResults pr = leo::compute(_appId, dataFile,
                                                 &_values, perColumnNormalization);
        _predictions = arma::conv_to<std::vector<double> >::from(pr.predictions);
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
          const Configuration& configuration,
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
            value = getMaximumBandwidth();
        }break;
        case PREDICTION_POWER:{
            value = getCurrentPower();
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

