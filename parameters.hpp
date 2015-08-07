/*
 * parameters.hpp
 *
 * Created on: 23/03/2015
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


#ifndef ADAPTIVE_FASTFLOW_PARAMETERS_HPP_
#define ADAPTIVE_FASTFLOW_PARAMETERS_HPP_

#include "external/rapidXml/rapidxml.hpp"

#include <cmath>
#include <fstream>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

namespace adpff{

class Observer;

using namespace std;
using namespace mammut::cpufreq;
using namespace mammut::topology;
using namespace mammut::utils;
using mammut::Communicator;
using mammut::Mammut;
using mammut::utils::enumStrings;

/// Possible contracts requested by the user.
typedef enum{
    // No contract required.
    CONTRACT_NONE = 0,

    // A specific utilization (expressed by two bounds [X, Y]) is requested.
    // Under this constraint, the configuration with minimum power consumption
    // is chosen.
    CONTRACT_PERF_UTILIZATION,

    // A specific minimum bandwidth is requested. Under this constraint, the
    // configuration with minimum power consumption is chosen.
    CONTRACT_PERF_BANDWIDTH,

    // A specific maximum completion time is requested. Under this constraint,
    // the configuration with minimum energy consumption is chosen.
    CONTRACT_PERF_COMPLETION_TIME,

    // A specific maximum cores power is requested. Under this constraint, the
    // configuration with best performance is chosen.
    CONTRACT_POWER_BUDGET
}ContractType;

template<> char const* enumStrings<ContractType>::data[] = {
    "NONE",
    "PERF_UTILIZATION",
    "PERF_BANDWIDTH",
    "PERF_COMPLETION_TIME",
    "POWER_BUDGET"
};


/// Possible reconfiguration strategies.
typedef enum{
    // Chooses YES or NO according to the possibility of the system.
    STRATEGY_FREQUENCY_AUTO = 0,

    // Reconfigures frequencies and number of workers. If contractType is
    // of type CONTRACT_PERF_*, it tries to minimize the consumed power.
    // If contractType is of type CONTRACT_POWER_*, it tries to maximize
    // the performances.
    STRATEGY_FREQUENCY_YES,

    // Does not reconfigure frequencies.
    STRATEGY_FREQUENCY_NO,

    // The frequencies are managed by OS governor.
    // The governor must be set in 'frequencyGovernor' parameter.
    STRATEGY_FREQUENCY_OS,

    // Reconfigures frequencies and number of workers. Tries always to
    // minimize the number of virtual cores used. Only valid if contract
    // type is CONTRACT_PERF_*.
    STRATEGY_FREQUENCY_MIN_CORES,
}StrategyFrequencies;

template<> char const* enumStrings<StrategyFrequencies>::data[] = {
    "YES",
    "NO",
    "OS",
    "MIN_CORES"
};

/// Possible mapping strategies.
typedef enum{
    // Mapping decisions will be performed by the operating system.
    STRATEGY_MAPPING_NO = 0,

    // Mapping strategy will be automatically chosen at runtime.
    STRATEGY_MAPPING_AUTO,

    // Tries to keep the threads as close as possible.
    STRATEGY_MAPPING_LINEAR,

    // Tries to make good use of the shared caches.
    // Particularly useful when threads have large working sets.
    STRATEGY_MAPPING_CACHE_EFFICIENT
}StrategyMapping;

template<> char const* enumStrings<StrategyMapping>::data[] = {
    "NO",
    "AUTO",
    "LINEAR",
    "CACHE_EFFICIENT"
};

/// Possible hyperthreading strategies.
typedef enum{
    // Hyperthreading is not used.
    STRATEGY_HT_NO = 0,

    // Hyperthreading is used since the beginning.
    STRATEGY_HT_YES_SOONER,

    // Hyperthreading is used only when we used all the physical cores.
    STRATEGY_HT_YES_LATER
}StrategyHyperthreading;

template<> char const* enumStrings<StrategyHyperthreading>::data[] = {
    "NO",
    "YES_SOONER",
    "YES_LATER"
};

/// Possible strategies to apply for unused virtual cores. For unused virtual
/// cores we mean those never used or those used only on some conditions.
typedef enum{
    // Nothing is done on unused virtual cores.
    STRATEGY_UNUSED_VC_NONE = 0,

    // Automatically choose one of the other strategies.
    STRATEGY_UNUSED_VC_AUTO,

    // Set the virtual cores to the lowest frequency (only
    // possible if all the other virtual cores on the same
    // domain are unused).
    STRATEGY_UNUSED_VC_LOWEST_FREQUENCY,

    // Turn off the virtual cores. They will not be anymore seen by the
    // operating system and it will not schedule anything on them.
    STRATEGY_UNUSED_VC_OFF
}StrategyUnusedVirtualCores;

template<> char const* enumStrings<StrategyUnusedVirtualCores>::data[] = {
    "NONE",
    "AUTO",
    "LOWEST_FREQUENCY",
    "OFF"
};

/// Possible strategies to use to predict power and performance values.
typedef enum{
    // Applies multivariate linear regression.
    STRATEGY_PREDICTION_REGRESSION_LINEAR = 0,

    // Applies a simple analytical model.
    STRATEGY_PREDICTION_SIMPLE
}StrategyPrediction;

template<> char const* enumStrings<StrategyPrediction>::data[] = {
    "SIMPLE",
    "REGRESSION_LINEAR"
};

/// Possible ways to detect if the model has an high error.
typedef enum{
    // Constant prediction error, set by the user.
    STRATEGY_PREDICTION_ERROR_CONSTANT = 0,

    // Prediction error set equal to the maximum between
    // coefficient of variation and maximum specified error.
    STRATEGY_PREDICTION_ERROR_COEFFVAR
}StrategyPredictionError;

template<> char const* enumStrings<StrategyPredictionError>::data[] = {
    "CONSTANT",
    "COEFFVAR"
};

/// Possible mappings for a service node (emitter or collector).
typedef enum{
    // The service node is mapped on a physical core where no workers
    // are mapped.
    SERVICE_NODE_MAPPING_ALONE = 0,

    // The service node is mapped on a physical core together with a worker.
    SERVICE_NODE_MAPPING_COLLAPSED,

    // The service node is mapped on an independent voltage domain and
    // is kept running at maximum performances (only available when
    // strategyFrequencies != STRATEGY_FREQUENCY_NO).
    SERVICE_NODE_MAPPING_PERFORMANCE
}ServiceNodeMapping;

template<> char const* enumStrings<ServiceNodeMapping>::data[] = {
    "ALONE",
    "COLLAPSED",
    "PERFORMANCE"
};

/// Possible ways to smooth the values.
typedef enum{
    // Simple moving average
    STRATEGY_SMOOTHING_MOVING_AVERAGE = 0,

    // Exponential moving average
    STRATEGY_SMOOTHING_EXPONENTIAL
}StrategySmoothing;

template<> char const* enumStrings<StrategySmoothing>::data[] = {
    "MOVING_AVERAGE",
    "EXPONENTIAL"
};

/// Possible ways to change the smoothing factor at runtime.
typedef enum{
    // Constant smoothing factor
    STRATEGY_SMOOTHING_FACTOR_CONST = 0,

    // Changes the smoothing factor according to the variation
    // of the data.
    STRATEGY_SMOOTHING_FACTOR_DYNAMIC
}StrategySmoothingFactor;

template<> char const* enumStrings<StrategySmoothingFactor>::data[] = {
    "CONST",
    "DYNAMIC"
};

/// Possible ways to select the calibration points.
typedef enum{
    // Random choice of the points.
    STRATEGY_CALIBRATION_RANDOM = 0,
    // Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992)
    STRATEGY_CALIBRATION_NIEDERREITER,
    // Antonov, Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980)
    STRATEGY_CALIBRATION_SOBOL,
    // J.H. Halton, Numerische Mathematik 2, 84-90 (1960) and B. Vandewoestyne
    // R. Cools Computational and Applied Mathematics 189, 1&2, 341-361 (2006)
    STRATEGY_CALIBRATION_HALTON,
    STRATEGY_CALIBRATION_HALTON_REVERSE,
}StrategyCalibration;

template<> char const* enumStrings<StrategyCalibration>::data[] = {
    "RANDOM",
    "NIEDERREITER",
    "SOBOL",
    "HALTON",
    "HALTON_REVERSE"
};

/// Possible ways to act when the manager finds no response on a queue
/// from a farm node.
typedef enum{
    // Spins.
    STRATEGY_POLLING_SPINNING = 0,

    // Suspend the execution on the processor for a short period of time.
    // However, the OS will see the thread as doing useful work and probably
    // will not get descheduled.
    STRATEGY_POLLING_PAUSE,

    // Sleeps for 0 nanoseconds. The OS may deschedule the thread.
    STRATEGY_POLLING_SLEEP_SMALL,

    // Sleeps for a period of time equal to the node average latency. The
    // OS may deschedule the thread.
    STRATEGY_POLLING_SLEEP_LATENCY
}StrategyPolling;

template<> char const* enumStrings<StrategyPolling>::data[] = {
    "SPINNING",
    "PAUSE",
    "SLEEP_SMALL",
    "SLEEP_LATENCY"
};

// How to decide if we have to consider new configurations
typedef enum{
    // Number of samples.
    STRATEGY_PERSISTENCE_SAMPLES = 0,

    // Number of tasks.
    STRATEGY_PERSISTENCE_TASKS,

    // Coefficient of variation.
    STRATEGY_PERSISTENCE_VARIATION
}StrategyPersistence;

template<> char const* enumStrings<StrategyPersistence>::data[] = {
    "SAMPLES",
    "TASKS",
    "VARIATION"
};

/// Possible parameters validation results.
typedef enum{
    // Parameters are ok.
    VALIDATION_OK = 0,

    // strategyFrequencies can be different from STRATEGY_FREQUENCY_NO
    // only if strategyMapping is different from STRATEGY_MAPPING_NO.
    VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING,

    // The specified frequency strategy is not supported on this machine.
    VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED,

    // Specified governor not supported on this machine.
    VALIDATION_GOVERNOR_UNSUPPORTED,

    // sensitiveEmitter or sensitiveCollector specified but frequency
    // strategy is STRATEGY_FREQUENCY_NO.
    VALIDATION_EC_SENSITIVE_WRONG_F_STRATEGY,

    // sensitiveEmitter or sensitiveCollector specified but highest
    // frequency can't be set.
    VALIDATION_EC_SENSITIVE_MISSING_GOVERNORS,

    // The bounds are invalid or the frequency strategy is not
    // STRATEGY_FREQUENCY_OS.
    VALIDATION_INVALID_FREQUENCY_BOUNDS,

    // Strategy for unused virtual cores requires turning off the virtual
    // cores but they can't be turned off.
    VALIDATION_UNUSED_VC_NO_OFF,

    // Strategy for unused virtual cores requires lowering the frequency but
    // frequency scaling not available.
    VALIDATION_UNUSED_VC_NO_FREQUENCIES,

    // Specified parameters are not valid for the specified contract.
    VALIDATION_WRONG_CONTRACT_PARAMETERS,

    // strategyFrequencies is STRATEGY_FREQUENCY_POWER_CONSERVATIVE but
    // the voltage file has not been specified or it does not exist.
    VALIDATION_VOLTAGE_FILE_NEEDED,

    // Fast reconfiguration not available.
    VALIDATION_NO_FAST_RECONF,
}AdaptivityParametersValidation;


class XmlContent{
private:
    char* _fileContentChars;
    rapidxml::xml_node<>* _root;
public:
    XmlContent(const string& fileName,
               const string& rootName){
        rapidxml::xml_document<> xmlContent;
        ifstream file(fileName.c_str());
        if(!file.is_open()){
            throw runtime_error("Impossible to read xml file " + fileName);
        }

        string fileContent;
        file.seekg(0, ios::end);
        fileContent.reserve(file.tellg());
        file.seekg(0, ios::beg);
        fileContent.assign((istreambuf_iterator<char>(file)),
                            istreambuf_iterator<char>());
        _fileContentChars = new char[fileContent.size() + 1];
        copy(fileContent.begin(), fileContent.end(), _fileContentChars);
        _fileContentChars[fileContent.size()] = '\0';
        xmlContent.parse<0>(_fileContentChars);
        _root = xmlContent.first_node(rootName.c_str());
    }

    ~XmlContent(){
        delete[] _fileContentChars;
    }

    void setBoolValue(const char* valueName, bool& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = string(node->value()).compare("true")?false:true;
        }
    }

    void setIntValue(const char* valueName, int& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = stringToInt(node->value());
        }
    }

    void setUintValue(const char* valueName, uint& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = stringToUint(node->value());
        }
    }

    void setUlongValue(const char* valueName, ulong& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = stringToUlong(node->value());
        }
    }

    void setDoubleValue(const char* valueName, double& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = stringToDouble(node->value());
        }
    }

    bool getStringValue(const char* valueName, string& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = node->value();
            return true;
        }else{
            return false;
        }
    }

    void setStringValue(const char* valueName, string& value){
        getStringValue(valueName, value);
    }

    template<typename T>
    void setEnumValue(const char* valueName, T& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            stringstream line(node->value());
            line >> enumFromStringInternal(value);
        }
    }
};

#define SETVALUE(XML, TYPE, NAME) XML.set##TYPE##Value(#NAME, NAME)

/*!
 * \class ArchData
 * \brief This class contains data specific to the architecture,
 *        obtained by microbenchmarks.
 *
 * This class contains data specific to the architecture,
 * obtained by microbenchmarks.
 */
typedef struct ArchData{
    // Number of ticks for a nanosecond.
    double ticksPerNs;

    // Number of ticks spent by a worker to reply to a monitoring request.
    double monitoringCost;

    // The file containing the voltage table. It is mandatory when
    // strategyFrequencies is STRATEGY_FREQUENCY_YES.
    string voltageTableFile;

    ArchData():ticksPerNs(0),
               monitoringCost(0),
               voltageTableFile(""){;}

    void loadXml(const string& archFileName){
        XmlContent xc(archFileName, "archData");
        SETVALUE(xc, Double, ticksPerNs);
        SETVALUE(xc, Double, monitoringCost);
        SETVALUE(xc, String, voltageTableFile);
    }
}ArchData;

/*!
 * \class AdaptivityParameters
 * \brief This class contains parameters for adaptivity choices.
 *
 * This class contains parameters for adaptivity choices.
 */
class AdaptivityParameters{

private:
    friend class ManagerFarm;

    /**
     * Sets default parameters
     */
    void setDefault(){
        contractType = CONTRACT_NONE;
        strategyMapping = STRATEGY_MAPPING_LINEAR;
        strategyHyperthreading = STRATEGY_HT_NO;
        strategyFrequencies = STRATEGY_FREQUENCY_AUTO;
        strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_NONE;
        strategyInactiveVirtualCores = STRATEGY_UNUSED_VC_NONE;
        strategyPrediction = STRATEGY_PREDICTION_REGRESSION_LINEAR;
        strategyPredictionErrorPrimary = STRATEGY_PREDICTION_ERROR_CONSTANT;
        strategyPredictionErrorSecondary = STRATEGY_PREDICTION_ERROR_CONSTANT;
        strategySmoothing = STRATEGY_SMOOTHING_EXPONENTIAL;
        strategyCalibration = STRATEGY_CALIBRATION_SOBOL;
        strategyPolling = STRATEGY_POLLING_SLEEP_LATENCY;
        strategyPersistence = STRATEGY_PERSISTENCE_SAMPLES;
        mappingEmitter = SERVICE_NODE_MAPPING_ALONE;
        mappingCollector = SERVICE_NODE_MAPPING_ALONE;
        frequencyGovernor = GOVERNOR_USERSPACE;
        turboBoost = false;
        frequencyLowerBound = 0;
        frequencyUpperBound = 0;
        fastReconfiguration = false;
        migrateCollector = false;
        smoothingFactor = 0;
        persistenceValue = 0;
        samplingInterval = 0;
        underloadThresholdFarm = 80.0;
        overloadThresholdFarm = 90.0;
        underloadThresholdWorker = 80.0;
        overloadThresholdWorker = 90.0;
        requiredBandwidth = 0;
        requiredCompletionTime = 0;
        expectedTasksNumber = 0;
        powerBudget = 0;
        maxPrimaryPredictionError = 5.0;
        maxSecondaryPredictionError = 5.0;
        maxMonitoringOverhead = 1.0;
        observer = NULL;
    }

    /**
     * Sets the default values for parameters that depends
     * from others.
     */
    void setDefaultPost(){
        if(!samplingInterval){
            //TODO Se questo sampling interval è molto minore della latenza
            // media di un task potrei settare il sampling interval alla
            // latenza media di un task.
            double msMonitoringCost = (archData.monitoringCost/
                                       archData.ticksPerNs*
                                       0.000001);
            samplingInterval = ceil(msMonitoringCost*
                                        (100.0 - maxMonitoringOverhead));
        }

        if(!smoothingFactor){
            switch(strategySmoothing){
                case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
                    smoothingFactor = 10;
                }break;
                case STRATEGY_SMOOTHING_EXPONENTIAL:{
                    smoothingFactor = 0.5;
                }break;
            }
        }

        if(!persistenceValue){
            switch(strategyPersistence){
                case STRATEGY_PERSISTENCE_SAMPLES:{
                    persistenceValue = 10;
                }break;
                case STRATEGY_PERSISTENCE_TASKS:{
                    persistenceValue = 1000;
                }break;
                case STRATEGY_PERSISTENCE_VARIATION:{
                    persistenceValue = 5;
                }break;
            }
        }
    }

    bool isGovernorAvailable(Governor g){
        return mammut.getInstanceCpuFreq()->isGovernorAvailable(g);
    }

    vector<Frequency> getAvailableFrequencies(){
        vector<Frequency> frequencies;
        CpuFreq* cpuFreq = mammut.getInstanceCpuFreq();
        if(cpuFreq){
            vector<Domain*> fDomains = cpuFreq->getDomains();
            if(fDomains.size()){
                frequencies = fDomains.front()->getAvailableFrequencies();
            }
        }
        return frequencies;
    }

    bool isUnusedVcOffAvailable(){
        vector<VirtualCore*> vc = mammut.getInstanceTopology()->
                                  getVirtualCores();
        for(size_t i = 0; i < vc.size(); i++){
            if(vc.at(i)->isHotPluggable()){
                return true;
            }
        }
    }

    bool isFrequencySettable(){
        vector<Frequency> frequencies = getAvailableFrequencies();
        return  isGovernorAvailable(GOVERNOR_USERSPACE) && frequencies.size();

    }

    bool isLowestFrequencySettable(){
        return isGovernorAvailable(GOVERNOR_POWERSAVE) ||
               isFrequencySettable();
    }

    bool isHighestFrequencySettable(){
        return isGovernorAvailable(GOVERNOR_PERFORMANCE) ||
               isFrequencySettable();
    }

    AdaptivityParametersValidation validateUnusedVc
                                   (StrategyUnusedVirtualCores& s){
        switch(s){
            case STRATEGY_UNUSED_VC_OFF:{
                if(!isUnusedVcOffAvailable()){
                    return VALIDATION_UNUSED_VC_NO_OFF;
                }
            }break;
            case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                if(!isLowestFrequencySettable()){
                    return VALIDATION_UNUSED_VC_NO_FREQUENCIES;
                }
            }break;
            case STRATEGY_UNUSED_VC_AUTO:{
                if(isUnusedVcOffAvailable()){
                    s = STRATEGY_UNUSED_VC_OFF;
                }else if(isLowestFrequencySettable()){
                    s = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY;
                }
            }
            default:
                break;
        }
        return VALIDATION_OK;
    }

    bool serviceNodePerformance(){
        return mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE ||
               mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE;
    }

    AdaptivityParametersValidation validateFrequencies(){
        vector<Frequency> availableFrequencies = getAvailableFrequencies();

        if(strategyFrequencies == STRATEGY_FREQUENCY_AUTO){
            if(isGovernorAvailable(GOVERNOR_USERSPACE) &&
               availableFrequencies.size()){
                strategyFrequencies = STRATEGY_FREQUENCY_YES;
            }else{
                strategyFrequencies = STRATEGY_FREQUENCY_NO;
            }
        }

        if(strategyFrequencies != STRATEGY_FREQUENCY_NO &&
           strategyMapping == STRATEGY_MAPPING_NO){
                return VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING;
        }

        if(strategyFrequencies == STRATEGY_FREQUENCY_YES){
            if(archData.voltageTableFile.empty() ||
               !existsFile(archData.voltageTableFile)){
                return VALIDATION_VOLTAGE_FILE_NEEDED;
            }
        }

        switch(strategyFrequencies){
            case STRATEGY_FREQUENCY_AUTO:{
                throw runtime_error("This should never happen.");
            }break;
            case STRATEGY_FREQUENCY_NO:{
                if(serviceNodePerformance()){
                    return VALIDATION_EC_SENSITIVE_WRONG_F_STRATEGY;
                }
            }break;
            case STRATEGY_FREQUENCY_OS:{
                if(!isGovernorAvailable(frequencyGovernor)){
                    return VALIDATION_GOVERNOR_UNSUPPORTED;
                }
                if(frequencyLowerBound || frequencyUpperBound){
                    if(!availableFrequencies.size()){
                        return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                    }

                    if(frequencyLowerBound){
                        if(!contains(availableFrequencies,
                                     frequencyLowerBound)){
                            return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                        }
                    }else{
                        frequencyLowerBound = availableFrequencies.front();
                    }

                    if(frequencyUpperBound){
                        if(!contains(availableFrequencies,
                                     frequencyUpperBound)){
                            return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                        }
                    }else{
                        frequencyUpperBound = availableFrequencies.back();
                    }
                }
                //TODO: Permettere di specificare i bound anche quando c'è
                //      FREQUENCY_YES
            }break;
            case STRATEGY_FREQUENCY_MIN_CORES:
            case STRATEGY_FREQUENCY_YES:{
                if(!availableFrequencies.size()){
                    return VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED;
                }

                frequencyGovernor = GOVERNOR_USERSPACE;
                if(!isGovernorAvailable(frequencyGovernor)){
                    return VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED;
                }

                if(serviceNodePerformance() &&
                   !isHighestFrequencySettable()){
                    return VALIDATION_EC_SENSITIVE_MISSING_GOVERNORS;
                }

            }break;
        }

        if(fastReconfiguration && !isHighestFrequencySettable()){
            return VALIDATION_NO_FAST_RECONF;
        }
        return VALIDATION_OK;
    }

    AdaptivityParametersValidation validateContract(){
        switch(contractType){
            case CONTRACT_PERF_UTILIZATION:{
                if(underloadThresholdFarm > overloadThresholdFarm     ||
                   underloadThresholdWorker > overloadThresholdWorker ||
                   underloadThresholdFarm < 0                         ||
                   overloadThresholdFarm > 100                        ||
                   underloadThresholdWorker < 0                       ||
                   overloadThresholdWorker > 100){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            case CONTRACT_PERF_BANDWIDTH:{
                if(requiredBandwidth <= 0){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            case CONTRACT_PERF_COMPLETION_TIME:{
                if(!expectedTasksNumber || !requiredCompletionTime){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            case CONTRACT_POWER_BUDGET:{
                if(powerBudget <= 0 ||
                   strategyFrequencies == STRATEGY_FREQUENCY_MIN_CORES ||
                   strategyPrediction == STRATEGY_PREDICTION_SIMPLE){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            default:{
                ;
            }break;
        }

        if(maxPrimaryPredictionError < 0      ||
           maxPrimaryPredictionError > 100.0  ||
           maxSecondaryPredictionError < 0    ||
           maxSecondaryPredictionError > 100.0){
            return VALIDATION_WRONG_CONTRACT_PARAMETERS;
        }

        return VALIDATION_OK;
    }

    void loadXml(const string& paramFileName){
        XmlContent xc(paramFileName, "adaptivityParameters");

        SETVALUE(xc, Enum, contractType);
        SETVALUE(xc, Enum, strategyMapping);
        SETVALUE(xc, Enum, strategyHyperthreading);
        SETVALUE(xc, Enum, strategyFrequencies);
        SETVALUE(xc, Enum, strategyUnusedVirtualCores);
        SETVALUE(xc, Enum, strategyInactiveVirtualCores);
        SETVALUE(xc, Enum, strategyPrediction);
        SETVALUE(xc, Enum, strategyPredictionErrorPrimary);
        SETVALUE(xc, Enum, strategyPredictionErrorSecondary);
        SETVALUE(xc, Enum, strategySmoothing);
        SETVALUE(xc, Enum, strategyCalibration);
        SETVALUE(xc, Enum, strategyPolling);
        SETVALUE(xc, Enum, strategyPersistence);
        SETVALUE(xc, Enum, mappingEmitter);
        SETVALUE(xc, Enum, mappingCollector);

        string g;
        if(xc.getStringValue("frequencyGovernor", g)){
            transform(g.begin(), g.end(), g.begin(), ::tolower);
            CpuFreq* cf = mammut.getInstanceCpuFreq();
            frequencyGovernor = cf->getGovernorFromGovernorName(g);
        }

        SETVALUE(xc, Bool, turboBoost);
        SETVALUE(xc, Uint, frequencyLowerBound);
        SETVALUE(xc, Uint, frequencyUpperBound);
        SETVALUE(xc, Bool, fastReconfiguration);
        SETVALUE(xc, Double, smoothingFactor);
        SETVALUE(xc, Double, persistenceValue);
        SETVALUE(xc, Uint, samplingInterval);
        SETVALUE(xc, Double, underloadThresholdFarm);
        SETVALUE(xc, Double, overloadThresholdFarm);
        SETVALUE(xc, Double, underloadThresholdWorker);
        SETVALUE(xc, Double, overloadThresholdWorker);
        SETVALUE(xc, Bool, migrateCollector);
        SETVALUE(xc, Double, requiredBandwidth);
        SETVALUE(xc, Uint, requiredCompletionTime);
        SETVALUE(xc, Ulong, expectedTasksNumber);
        SETVALUE(xc, Double, powerBudget);
        SETVALUE(xc, Double, maxPrimaryPredictionError);
        SETVALUE(xc, Double, maxSecondaryPredictionError);
        SETVALUE(xc, Double, maxMonitoringOverhead);
    }

public:
    // The mammut modules handler.
    Mammut mammut;

    // Architecture's specific data.
    ArchData archData;

    // The contract type that must be respected by the application
    // [default = CONTRACT_NONE].
    ContractType contractType;

    //  The mapping strategy [default = STRATEGY_MAPPING_LINEAR].
    StrategyMapping strategyMapping;

    // The hyperthreading strategy [default = STRATEGY_HT_NO].
    StrategyHyperthreading strategyHyperthreading;

    // The frequency strategy. It can be different from STRATEGY_FREQUENCY_NO
    // only if strategyMapping is different from STRATEGY_MAPPING_NO
    // [default = STRATEGY_FREQUENCY_AUTO].
    StrategyFrequencies strategyFrequencies;

    // Strategy for virtual cores that are never used
    // [default = STRATEGY_UNUSED_VC_NONE].
    StrategyUnusedVirtualCores strategyUnusedVirtualCores;

    // Strategy for virtual cores that become inactive after a workers
    // reconfiguration [default = STRATEGY_UNUSED_VC_NONE].
    StrategyUnusedVirtualCores strategyInactiveVirtualCores;

    // Strategy to be used to predict power and performance values
    // [default = STRATEGY_PREDICTION_REGRESSION_LINEAR].
    StrategyPrediction strategyPrediction;

    // Strategy to be used for toleration of prediction
    // errors on primary value [default = STRATEGY_PREDICTION_ERROR_CONSTANT].
    StrategyPredictionError strategyPredictionErrorPrimary;

    // Strategy to be used for toleration of prediction
    // errors on secondary value [default = STRATEGY_PREDICTION_ERROR_CONSTANT].
    StrategyPredictionError strategyPredictionErrorSecondary;

    // Smoothing strategy [default = STRATEGY_SMOOTHING_EXPONENTIAL].
    StrategySmoothing strategySmoothing;

    // Calibration strategy [default = STRATEGY_CALIBRATION_SOBOL].
    StrategyCalibration strategyCalibration;

    // Polling strategy [default = STRATEGY_POLLING_SLEEP_LATENCY].
    StrategyPolling strategyPolling;

    // Persistence strategy [default = STRATEGY_PERSISTENCE_SAMPLES].
    StrategyPersistence strategyPersistence;

    // The frequency governor (only used when strategyFrequencies is
    // STRATEGY_FREQUENCY_OS) [default = GOVERNOR_USERSPACE].
    Governor frequencyGovernor;

    // Flag to enable/disable cores turbo boosting [default = false].
    bool turboBoost;

    // The frequency lower bound (only if strategyFrequency is
    // STRATEGY_FREQUENCY_OS) [default = unused].
    Frequency frequencyLowerBound;

    // The frequency upper bound (only if strategyFrequency is
    // STRATEGY_FREQUENCY_OS) [default = unused].
    Frequency frequencyUpperBound;

    // If true, before changing the number of workers the frequency will be
    // set to maximum to reduce the latency of the reconfiguration. The
    // frequency will be be set again to the correct value after the farm
    // is restarted [default = false].
    bool fastReconfiguration;

    // Emitter mapping [default = SERVICE_NODE_MAPPING_ALONE].
    ServiceNodeMapping mappingEmitter;

    // Collector mapping [default = SERVICE_NODE_MAPPING_ALONE].
    ServiceNodeMapping mappingCollector;

    // If true, when a reconfiguration occur, the collector is migrated to a
    // different virtual core (if needed) [default = false].
    bool migrateCollector;

    // The smoothing factor. It's meaning changes according to the smoothing
    // strategy adopted. If 0, default value will be set.
    // [default = 10 for moving average, 0.5 for exponential].
    double smoothingFactor;

    // The persistence value. It's meaning changes according to the persistence
    // strategy adopted. If 0, default value will be set.
    // [default = 10 for samples, 1000 for tasks, 5 for variation].
    double persistenceValue;

    // The length of the sampling interval (in milliseconds) for the data
    // reading. If 0, it will be automatically computed such to have a low
    // overhead [default = 0].
    uint32_t samplingInterval;

    // The underload threshold for the entire farm. It is valid only if
    // contractType is CONTRACT_UTILIZATION [default = 80.0].
    double underloadThresholdFarm;

    // The overload threshold for the entire farm. It is valid only if
    // contractType is CONTRACT_UTILIZATION [default = 90.0].
    double overloadThresholdFarm;

    // The underload threshold for a single worker. It is valid only if
    // contractType is CONTRACT_UTILIZATION [default = 80.0].
    double underloadThresholdWorker;

    // The overload threshold for a single worker. It is valid only if
    // contractType is CONTRACT_UTILIZATION [default = 90.0].
    double overloadThresholdWorker;

    // The bandwidth required for the application (expressed as tasks/sec).
    // It is valid only if contractType is CONTRACT_BANDWIDTH
    // [default = unused].
    double requiredBandwidth;

    // The required completion time for the application (in seconds). It is
    // valid only if contractType is CONTRACT_COMPLETION_TIME
    // [default = unused].
    uint requiredCompletionTime;

    // The number of task expected for this computation. It is
    // valid only if contractType is CONTRACT_COMPLETION_TIME
    // [default = unused].
    ulong expectedTasksNumber;

    // The maximum cores power to be used. It is
    // valid only if contractType is CONTRACT_POWER_BUDGET [default = unused].
    double powerBudget;

    // Maximum error percentage allowed for prediction of primary
    // value. If 0, then it will be set equal to the coefficient
    // of variation of the primary value [default = 5.0].
    double maxPrimaryPredictionError;

    // Maximum error percentage allowed for prediction of secondary
    // value. If 0, then it will be set equal to the coefficient
    // of variation of the secondary value [default = 5.0].
    double maxSecondaryPredictionError;

    // The maximum percentage of monitoring overhead, in the range (0, 100).
    // [default = 1.0].
    double maxMonitoringOverhead;

    // The observer object. It will be called every samplingInterval
    // milliseconds to monitor the adaptivity behaviour [default = NULL].
    Observer* observer;

    /**
     * Creates the adaptivity parameters.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    AdaptivityParameters(Communicator* const communicator = NULL):
          mammut(communicator){
        setDefault();
    }

    /**
     * Creates the adaptivity parameters.
     * @param paramFileName The name of the XML file containing the adaptivity
     *        parameters.
     * @param archFileName The name of the XML file containing the
     *        architectural data.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    AdaptivityParameters(const string& paramFileName,
                         const string& archFileName,
                         Communicator* const communicator = NULL):
          mammut(communicator){
        setDefault();
        archData.loadXml(archFileName);
        loadXml(paramFileName);
    }


    /**
     * Destroyes these parameters.
     */
    ~AdaptivityParameters(){
        ;
    }

    /**
     * Validates these parameters.
     * @return The validation result.
     */
    AdaptivityParametersValidation validate(){
        AdaptivityParametersValidation r = VALIDATION_OK;
        setDefaultPost();

        vector<VirtualCore*> virtualCores;
        virtualCores = mammut.getInstanceTopology()->getVirtualCores();

        /** Validate frequency strategies. **/
        r = validateFrequencies();
        if(r != VALIDATION_OK){return r;}

        /** Validate unused cores strategy. **/
        r = validateUnusedVc(strategyInactiveVirtualCores);
        if(r != VALIDATION_OK){return r;}
        validateUnusedVc(strategyUnusedVirtualCores);
        if(r != VALIDATION_OK){return r;}

        /** Validate contract parameters. **/
        r = validateContract();
        if(r != VALIDATION_OK){return r;}

        /** Validate unsupported strategies. **/
        if(strategyHyperthreading == STRATEGY_HT_YES_SOONER ||
           strategyMapping == STRATEGY_MAPPING_CACHE_EFFICIENT){
            throw runtime_error("Not yet supported.");
        }

        return VALIDATION_OK;
    }
};

}

#endif /* ADAPTIVE_FASTFLOW_PARAMETERS_HPP_ */
