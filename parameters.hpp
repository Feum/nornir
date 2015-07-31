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

#include <fstream>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

using mammut::utils::enumStrings;

namespace adpff{

class Observer;
using namespace mammut;

/// Possible contracts requested by the user.
typedef enum{
    CONTRACT_NONE = 0, ///< No contract required.
    CONTRACT_PERF_UTILIZATION, ///< A specific utilization (expressed by two bounds [X, Y]) is requested.
                               ///< Under this constraint, the configuration with minimum power consumption
                               ///< is chosen.
    CONTRACT_PERF_BANDWIDTH, ///< A specific minimum bandwidth is requested. Under this constraint, the
                             ///< configuration with minimum power consumption is chosen.
    CONTRACT_PERF_COMPLETION_TIME, ///< A specific maximum completion time is requested. Under this constraint, the
                                   ///< configuration with minimum energy consumption is chosen.
    CONTRACT_POWER_BUDGET ///< A specific maximum cores power is requested. Under this constraint, the configuration
                          ///< with best performance is chosen.
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
    STRATEGY_FREQUENCY_YES = 0, ///< Reconfigures frequencies and number of workers. If contractType is
                                ///< of type CONTRACT_PERF_*, it tries to minimize the consumed power.
                                ///< If contractType is of type CONTRACT_POWER_*, it tries to maximize
                                ///< the performances.
    STRATEGY_FREQUENCY_NO, ///< Does not reconfigure frequencies.
    STRATEGY_FREQUENCY_OS, ///< Reconfigures the number of workers. The frequencies are managed by OS governor.
                           ///< The governor must be specified with 'setFrequengyGovernor' call.
    STRATEGY_FREQUENCY_MIN_CORES, ///< Reconfigures frequencies and number of workers. Tries always to
                                  ///< minimize the number of virtual cores used. Only valid if contract
                                  ///< type is CONTRACT_PERF_*.

}StrategyFrequencies;

template<> char const* enumStrings<StrategyFrequencies>::data[] = {
    "YES",
    "NO",
    "OS",
    "MIN_CORES"
};

/// Possible mapping strategies.
typedef enum{
    STRATEGY_MAPPING_NO = 0, ///< Mapping decisions will be performed by the operating system.
    STRATEGY_MAPPING_AUTO, ///< Mapping strategy will be automatically chosen at runtime.
    STRATEGY_MAPPING_LINEAR, ///< Tries to keep the threads as close as possible.
    STRATEGY_MAPPING_CACHE_EFFICIENT ///< Tries to make good use of the shared caches.
                                     ///< Particularly useful when threads have large working sets.
}StrategyMapping;

template<> char const* enumStrings<StrategyMapping>::data[] = {
    "NO",
    "AUTO",
    "LINEAR",
    "CACHE_EFFICIENT"
};

/// Possible hyperthreading strategies.
typedef enum{
    STRATEGY_HT_NO = 0, ///< Hyperthreading is not used.
    STRATEGY_HT_YES_SOONER, ///< Hyperthreading is used since the beginning.
    STRATEGY_HT_YES_LATER ///< Hyperthreading is used only when we used all the physical cores.
}StrategyHyperthreading;

template<> char const* enumStrings<StrategyHyperthreading>::data[] = {
    "NO",
    "YES_SOONER",
    "YES_LATER"
};

/// Possible strategies to apply for unused virtual cores. For unused virtual cores we mean
/// those never used or those used only on some conditions.
typedef enum{
    STRATEGY_UNUSED_VC_NONE = 0, ///< Nothing is done on unused virtual cores.
    STRATEGY_UNUSED_VC_AUTO, ///< Automatically choose one of the other strategies.
    STRATEGY_UNUSED_VC_LOWEST_FREQUENCY, ///< Set the virtual cores to the lowest frequency (only
                                         ///< possible if all the other virtual cores on the same
                                         ///< domain are unused).
    STRATEGY_UNUSED_VC_OFF ///< Turn off the virtual cores. They will not be anymore seen by the
                           ///< operating system and it will not schedule anything on them.
}StrategyUnusedVirtualCores;

template<> char const* enumStrings<StrategyUnusedVirtualCores>::data[] = {
    "NONE",
    "AUTO",
    "LOWEST_FREQUENCY",
    "OFF"
};

/// Possible strategies to use to predict power and performance values.
typedef enum{
    STRATEGY_PREDICTION_SIMPLE = 0,       ///< Applies a simple analytical model.
    STRATEGY_PREDICTION_REGRESSION_LINEAR ///< Applies multivariate linear regression.
}StrategyPrediction;

template<> char const* enumStrings<StrategyPrediction>::data[] = {
    "SIMPLE",
    "REGRESSION_LINEAR"
};

/// Possible mappings for a service node (emitter or collector).
typedef enum{
    SERVICE_NODE_MAPPING_ALONE = 0, ///< The service node is mapped on a physical core where no workers are mapped.
    SERVICE_NODE_MAPPING_COLLAPSED, ///< The service node is mapped on a physical core together with a worker.
    SERVICE_NODE_MAPPING_PERFORMANCE ///< The service node is mapped on an independent voltage domain and is kept
                                     ///< running at maximum performances (only available when
                                     ///< strategyFrequencies != STRATEGY_FREQUENCY_NO).
}ServiceNodeMapping;

template<> char const* enumStrings<ServiceNodeMapping>::data[] = {
    "ALONE",
    "COLLAPSED",
    "PERFORMANCE"
};

/// Possible ways to smooth the values.
typedef enum{
    STRATEGY_SMOOTHING_MOVING_AVERAGE = 0, ///< Simple moving average
    STRATEGY_SMOOTHING_EXPONENTIAL ///< Exponential moving average
}StrategySmoothing;

template<> char const* enumStrings<StrategySmoothing>::data[] = {
    "MOVING_AVERAGE",
    "EXPONENTIAL"
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

/// Possible parameters validation results.
typedef enum{
    VALIDATION_OK = 0, ///< Parameters are ok.
    VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING, ///< strategyFrequencies can be different from STRATEGY_FREQUENCY_NO
                                                    ///< only if strategyMapping is different from STRATEGY_MAPPING_NO.
    VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED, ///< The specified frequency strategy is not supported
                                               ///< on this machine.
    VALIDATION_GOVERNOR_UNSUPPORTED, ///< Specified governor not supported on this machine.
    VALIDATION_STRATEGY_MAPPING_UNSUPPORTED, ///< Specified mapping strategy not supported on this machine.
    VALIDATION_EC_SENSITIVE_WRONG_F_STRATEGY, ///< sensitiveEmitter or sensitiveCollector specified but frequency
                                              ///< strategy is STRATEGY_FREQUENCY_NO.
    VALIDATION_EC_SENSITIVE_MISSING_GOVERNORS, ///< sensitiveEmitter or sensitiveCollector specified but highest
                                               ///< frequency can't be set.
    VALIDATION_INVALID_FREQUENCY_BOUNDS, ///< The bounds are invalid or the frequency strategy is not STRATEGY_FREQUENCY_OS.
    VALIDATION_UNUSED_VC_NO_OFF, ///< Strategy for unused virtual cores requires turning off the virtual cores but they
                                 ///< can't be turned off.
    VALIDATION_UNUSED_VC_NO_FREQUENCIES, ///< Strategy for unused virtual cores requires lowering the frequency but
                                         ///< frequency scaling not available.
    VALIDATION_WRONG_CONTRACT_PARAMETERS, ///< Specified parameters are not valid for the specified contract.
    VALIDATION_VOLTAGE_FILE_NEEDED, ///< strategyFrequencies is STRATEGY_FREQUENCY_POWER_CONSERVATIVE but the voltage file
                                    ///< has not been specified or it does not exist.
    VALIDATION_NO_FAST_RECONF, ///< Fast reconfiguration not available.
}AdaptivityParametersValidation;


class XmlContent{
private:
    char* _fileContentChars;
    rapidxml::xml_node<>* _root;
public:
    XmlContent(const std::string& fileName,
               const std::string& rootName){
        rapidxml::xml_document<> xmlContent;
        std::ifstream file(fileName.c_str());
        if(!file.is_open()){
            throw std::runtime_error("Impossible to read xml file " + fileName);
        }

        std::string fileContent;
        file.seekg(0, std::ios::end);
        fileContent.reserve(file.tellg());
        file.seekg(0, std::ios::beg);
        fileContent.assign((std::istreambuf_iterator<char>(file)),
                            std::istreambuf_iterator<char>());
        _fileContentChars = new char[fileContent.size() + 1];
        std::copy(fileContent.begin(), fileContent.end(), _fileContentChars);
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
            value = std::string(node->value()).compare("true")?false:true;
        }
    }

    void setIntValue(const char* valueName, int& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = utils::stringToInt(node->value());
        }
    }

    void setUintValue(const char* valueName, uint& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = utils::stringToUint(node->value());
        }
    }

    void setUlongValue(const char* valueName, ulong& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = utils::stringToUlong(node->value());
        }
    }

    void setDoubleValue(const char* valueName, double& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = utils::stringToDouble(node->value());
        }
    }

    bool getStringValue(const char* valueName, std::string& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            value = node->value();
            return true;
        }else{
            return false;
        }
    }

    void setStringValue(const char* valueName, std::string& value){
        getStringValue(valueName, value);
    }

    template<typename T>
    void setEnumValue(const char* valueName, T& value){
        rapidxml::xml_node<> *node = NULL;
        node = _root->first_node(valueName);
        if(node){
            std::stringstream line(node->value());
            line >> mammut::utils::enumFromStringInternal(value);
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

    // The file containing the voltage table. It is mandatory when
    // strategyFrequencies is STRATEGY_FREQUENCY_YES.
    std::string voltageTableFile;

    ArchData():ticksPerNs(0),
               voltageTableFile(""){;}

    void loadXml(const std::string& archFileName){
        XmlContent xc(archFileName, "archData");
        SETVALUE(xc, Double, ticksPerNs);
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
    template<typename lb_t, typename gt_t>
    friend class adp_ff_farm;

    friend class AdaptivityManagerFarm;

    /**
     * Sets default parameters
     */
    void setDefault(){
        contractType = CONTRACT_NONE;
        strategyMapping = STRATEGY_MAPPING_LINEAR;
        strategyHyperthreading = STRATEGY_HT_NO;
        strategyFrequencies = STRATEGY_FREQUENCY_NO;
        strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_NONE;
        strategyInactiveVirtualCores = STRATEGY_UNUSED_VC_NONE;
        strategyPrediction = STRATEGY_PREDICTION_SIMPLE;
        strategySmoothing = STRATEGY_SMOOTHING_MOVING_AVERAGE;
        strategyCalibration = STRATEGY_CALIBRATION_SOBOL;
        strategyPolling = STRATEGY_POLLING_PAUSE;
        mappingEmitter = SERVICE_NODE_MAPPING_ALONE;
        mappingCollector = SERVICE_NODE_MAPPING_ALONE;
        frequencyGovernor = cpufreq::GOVERNOR_USERSPACE;
        turboBoost = false;
        frequencyLowerBound = 0;
        frequencyUpperBound = 0;
        fastReconfiguration = false;
        migrateCollector = false;
        numSamples = 10;
        alphaExpAverage = 0.5;
        samplingInterval = 1000.0;
        underloadThresholdFarm = 80.0;
        overloadThresholdFarm = 90.0;
        underloadThresholdWorker = 80.0;
        overloadThresholdWorker = 90.0;
        requiredBandwidth = 0;
        requiredCompletionTime = 0;
        expectedTasksNumber = 0;
        powerBudget = 0;
        maxPredictionError = 10.0;
        observer = NULL;
    }

    void loadXml(const std::string& paramFileName){
        XmlContent xc(paramFileName, "adaptivityParameters");

        SETVALUE(xc, Enum, contractType);
        SETVALUE(xc, Enum, strategyMapping);
        SETVALUE(xc, Enum, strategyHyperthreading);
        SETVALUE(xc, Enum, strategyFrequencies);
        SETVALUE(xc, Enum, strategyUnusedVirtualCores);
        SETVALUE(xc, Enum, strategyInactiveVirtualCores);
        SETVALUE(xc, Enum, strategyPrediction);
        SETVALUE(xc, Enum, strategySmoothing);
        SETVALUE(xc, Enum, strategyCalibration);
        SETVALUE(xc, Enum, strategyPolling);
        SETVALUE(xc, Enum, mappingEmitter);
        SETVALUE(xc, Enum, mappingCollector);

        std::string g;
        if(xc.getStringValue("frequencyGovernor", g)){
            std::transform(g.begin(), g.end(), g.begin(), ::tolower);
            mammut::cpufreq::CpuFreq* cf = mammut.getInstanceCpuFreq();
            frequencyGovernor = cf->getGovernorFromGovernorName(g);
        }

        SETVALUE(xc, Bool, turboBoost);
        SETVALUE(xc, Uint, frequencyLowerBound);
        SETVALUE(xc, Uint, frequencyUpperBound);
        SETVALUE(xc, Bool, fastReconfiguration);
        SETVALUE(xc, Uint, numSamples);
        SETVALUE(xc, Double, alphaExpAverage);
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
        SETVALUE(xc, Double, maxPredictionError);
    }

public:
    // The mammut modules handler.
    mammut::Mammut mammut;

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
    // [default = STRATEGY_FREQUENCY_NO].
    StrategyFrequencies strategyFrequencies;

    // Strategy for virtual cores that are never used
    // [default = STRATEGY_UNUSED_VC_NONE].
    StrategyUnusedVirtualCores strategyUnusedVirtualCores;

    // Strategy for virtual cores that become inactive after a workers
    // reconfiguration [default = STRATEGY_UNUSED_VC_NONE].
    StrategyUnusedVirtualCores strategyInactiveVirtualCores;

    // Strategy to be used to predict power and performance values
    // [default = STRATEGY_PREDICTION_SIMPLE].
    StrategyPrediction strategyPrediction;

    // Smoothing strategy [default = STRATEGY_SMOOTHING_SIMPLE_AVERAGE].
    StrategySmoothing strategySmoothing;

    // Calibration strategy [default = STRATEGY_CALIBRATION_SOBOL].
    StrategyCalibration strategyCalibration;

    // Polling strategy [default = STRATEGY_POLLING_PAUSE].
    StrategyPolling strategyPolling;

    // The frequency governor (only used when strategyFrequencies is
    // STRATEGY_FREQUENCY_OS) [default = GOVERNOR_USERSPACE].
    cpufreq::Governor frequencyGovernor;

    // Flag to enable/disable cores turbo boosting [default = false].
    bool turboBoost;

    // The frequency lower bound (only if strategyFrequency is
    // STRATEGY_FREQUENCY_OS) [default = unused].
    cpufreq::Frequency frequencyLowerBound;

    // The frequency upper bound (only if strategyFrequency is
    // STRATEGY_FREQUENCY_OS) [default = unused].
    cpufreq::Frequency frequencyUpperBound;

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

    // The minimum number of samples used to take reconfiguration
    // decisions [default = 10].
    uint32_t numSamples;

    // The alpha to be used in exponential moving average [default = 0.5].
    double alphaExpAverage;

    // The length of the sampling interval (in milliseconds) for the data
    // reading. If 0, it will be automatically computed such to have a low
    // overhead [default = 0].
    uint32_t samplingInterval;

    double underloadThresholdFarm; ///< The underload threshold for the entire farm. It is valid only if
                                   ///< contractType is CONTRACT_UTILIZATION [default = 80.0].
    double overloadThresholdFarm; ///< The overload threshold for the entire farm. It is valid only if
                                  ///< contractType is CONTRACT_UTILIZATION [default = 90.0].
    double underloadThresholdWorker; ///< The underload threshold for a single worker. It is valid only if
                                     ///< contractType is CONTRACT_UTILIZATION [default = 80.0].
    double overloadThresholdWorker; ///< The overload threshold for a single worker. It is valid only if
                                    ///< contractType is CONTRACT_UTILIZATION [default = 90.0].
    double requiredBandwidth; ///< The bandwidth required for the application (expressed as tasks/sec).
                              ///< It is valid only if contractType is CONTRACT_BANDWIDTH [default = unused].
    uint requiredCompletionTime; ///< The required completion time for the application (in seconds). It is
                                 ///< valid only if contractType is CONTRACT_COMPLETION_TIME [default = unused].
    ulong expectedTasksNumber; ///< The number of task expected for this computation. It is
                               ///< valid only if contractType is CONTRACT_COMPLETION_TIME [default = unused].
    double powerBudget;           ///< The maximum cores power to be used. It is
                                  ///< valid only if contractType is CONTRACT_POWER_BUDGET [default = unused].
    double maxPredictionError; ///< Maximum error percentage allowed for prediction [default = 10.0].
    Observer* observer; ///< The observer object. It will be called every samplingInterval milliseconds
                                   ///< to monitor the adaptivity behaviour [default = NULL].

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
    AdaptivityParameters(const std::string& paramFileName,
                         const std::string& archFileName,
                         Communicator* const communicator = NULL):
          mammut(communicator){
        setDefault();
        loadXml(paramFileName);
        archData.loadXml(archFileName);
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
        std::vector<cpufreq::Domain*> frequencyDomains = mammut.getInstanceCpuFreq()->getDomains();
        std::vector<topology::VirtualCore*> virtualCores = mammut.getInstanceTopology()->getVirtualCores();
        std::vector<cpufreq::Frequency> availableFrequencies;
        if(frequencyDomains.size()){
            availableFrequencies = frequencyDomains.at(0)->getAvailableFrequencies();
        }

        if(strategyFrequencies != STRATEGY_FREQUENCY_NO && strategyMapping == STRATEGY_MAPPING_NO){
            return VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING;
        }

        /** Validate frequency strategies. **/
        if(strategyFrequencies != STRATEGY_FREQUENCY_NO){
            if(!frequencyDomains.size()){
                return VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED;
            }

            if(strategyFrequencies != STRATEGY_FREQUENCY_OS){
                frequencyGovernor = cpufreq::GOVERNOR_USERSPACE;
                if(!mammut.getInstanceCpuFreq()->isGovernorAvailable(frequencyGovernor)){
                    return VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED;
                }
            }
            if(((mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE) || (mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE)) &&
               !mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_PERFORMANCE) &&
               !mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_USERSPACE)){
                return VALIDATION_EC_SENSITIVE_MISSING_GOVERNORS;
            }
        }else{
            if((mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE) || (mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE)){
                return VALIDATION_EC_SENSITIVE_WRONG_F_STRATEGY;
            }
        }

        /** Validate governor availability. **/
        if(!mammut.getInstanceCpuFreq()->isGovernorAvailable(frequencyGovernor)){
            return VALIDATION_GOVERNOR_UNSUPPORTED;
        }

        /** Validate mapping strategy. **/
        if(strategyMapping == STRATEGY_MAPPING_CACHE_EFFICIENT){
            return VALIDATION_STRATEGY_MAPPING_UNSUPPORTED;
        }

        /** Validate frequency bounds. **/
        if(frequencyLowerBound || frequencyUpperBound){
            if(strategyFrequencies == STRATEGY_FREQUENCY_OS){
                if(!availableFrequencies.size()){
                    return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                }

                if(frequencyLowerBound){
                    if(!utils::contains(availableFrequencies, frequencyLowerBound)){
                        return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                    }
                }else{
                    frequencyLowerBound = availableFrequencies.at(0);
                }

                if(frequencyUpperBound){
                    if(!utils::contains(availableFrequencies, frequencyUpperBound)){
                        return VALIDATION_INVALID_FREQUENCY_BOUNDS;
                    }
                }else{
                    frequencyUpperBound = availableFrequencies.back();
                }
            }else{
                return VALIDATION_INVALID_FREQUENCY_BOUNDS;
            }
        }

        /** Validate unused cores strategy. **/
        switch(strategyInactiveVirtualCores){
            case STRATEGY_UNUSED_VC_OFF:{
                bool hotPluggableFound = false;
                for(size_t i = 0; i < virtualCores.size(); i++){
                    if(virtualCores.at(i)->isHotPluggable()){
                        hotPluggableFound = true;
                    }
                }
                if(!hotPluggableFound){
                    return VALIDATION_UNUSED_VC_NO_OFF;
                }
            }break;
            case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                if(!mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_POWERSAVE) &&
                   !mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_USERSPACE)){
                    return VALIDATION_UNUSED_VC_NO_FREQUENCIES;
                }
            }break;
            default:
                break;
        }

        /** Validate contract parameters. **/
        switch(contractType){
            case CONTRACT_PERF_UTILIZATION:{
                if((underloadThresholdFarm > overloadThresholdFarm) ||
                   (underloadThresholdWorker > overloadThresholdWorker) ||
                   underloadThresholdFarm < 0 || overloadThresholdFarm > 100 ||
                   underloadThresholdWorker < 0 || overloadThresholdWorker > 100){
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
                if(powerBudget <= 0 || (strategyFrequencies == STRATEGY_FREQUENCY_MIN_CORES) ||
                   strategyPrediction == STRATEGY_PREDICTION_SIMPLE){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            default:{
                ;
            }break;
        }

        if(maxPredictionError < 0 || maxPredictionError > 100.0){
            return VALIDATION_WRONG_CONTRACT_PARAMETERS;
        }

        /** Validate voltage table. **/
        if(strategyFrequencies == STRATEGY_FREQUENCY_YES){
            if(archData.voltageTableFile.empty() ||
               !utils::existsFile(archData.voltageTableFile)){
                return VALIDATION_VOLTAGE_FILE_NEEDED;
            }
        }

        /** Validate fast reconfiguration. **/
        if(fastReconfiguration){
            if(!mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_PERFORMANCE) &&
               (!mammut.getInstanceCpuFreq()->isGovernorAvailable(cpufreq::GOVERNOR_USERSPACE) ||
                !availableFrequencies.size())){
                return VALIDATION_NO_FAST_RECONF;
            }
        }

        /** Validate Hyperthreading strategies. **/
        if(strategyHyperthreading == STRATEGY_HT_YES_SOONER){
            throw std::runtime_error("Not yet supported.");
        }

        return VALIDATION_OK;
    }
};

}

#endif /* ADAPTIVE_FASTFLOW_PARAMETERS_HPP_ */
