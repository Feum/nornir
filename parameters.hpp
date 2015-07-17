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

namespace adpff{

class adp_ff_farm_observer;

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

/// Possible mapping strategies.
typedef enum{
    STRATEGY_MAPPING_NO = 0, ///< Mapping decisions will be performed by the operating system.
    STRATEGY_MAPPING_AUTO, ///< Mapping strategy will be automatically chosen at runtime.
    STRATEGY_MAPPING_LINEAR, ///< Tries to keep the threads as close as possible.
    STRATEGY_MAPPING_CACHE_EFFICIENT ///< Tries to make good use of the shared caches.
                                     ///< Particularly useful when threads have large working sets.
}StrategyMapping;

/// Possible hyperthreading strategies.
typedef enum{
    STRATEGY_HT_NO = 0, ///< Hyperthreading is not used.
    STRATEGY_HT_YES_SOONER, ///< Hyperthreading is used since the beginning.
    STRATEGY_HT_YES_LATER ///< Hyperthreading is used only when we used all the physical cores.
}StrategyHyperthreading;

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

/// Possible strategies to use to predict power and performance values.
typedef enum{
    STRATEGY_PREDICTION_SIMPLE = 0,       ///< Applies a simple analytical model.
    STRATEGY_PREDICTION_REGRESSION_LINEAR ///< Applies multivariate linear regression.
}StrategyPrediction;

/// Possible mappings for a service node (emitter or collector).
typedef enum{
    SERVICE_NODE_MAPPING_ALONE = 0, ///< The service node is mapped on a physical core where no workers are mapped.
    SERVICE_NODE_MAPPING_COLLAPSED, ///< The service node is mapped on a physical core together with a worker.
    SERVICE_NODE_MAPPING_PERFORMANCE ///< The service node is mapped on an independent voltage domain and is kept
                                     ///< running at maximum performances (only available when
                                     ///< strategyFrequencies != STRATEGY_FREQUENCY_NO).
}ServiceNodeMapping;

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
                                               ///< frequency can't be set..
    VALIDATION_INVALID_FREQUENCY_BOUNDS, ///< The bounds are invalid or the frequency strategy is not STRATEGY_FREQUENCY_OS.
    VALIDATION_UNUSED_VC_NO_OFF, ///< Strategy for unused virtual cores requires turning off the virtual cores but they
                                 ///< can't be turned off.
    VALIDATION_UNUSED_VC_NO_FREQUENCIES, ///< Strategy for unused virtual cores requires lowering the frequency but
                                         ///< frequency scaling not available.
    VALIDATION_WRONG_CONTRACT_PARAMETERS, ///< Specified parameters are not valid for the specified contract.
    VALIDATION_VOLTAGE_FILE_NEEDED, ///< strategyFrequencies is STRATEGY_FREQUENCY_POWER_CONSERVATIVE but the voltage file
                                    ///< has not been specified or it does not exist.
    VALIDATION_NO_FAST_RECONF, ///< Fast reconfiguration not available.
    VALIDATION_WRONG_REGRESSION_PARAMETERS ///< If strategyPrediction is STRATEGY_PREDICTION_REGRESSION_LINEAR but the
                                           ///< parameters are wrong.
}AdaptivityParametersValidation;

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
        mappingEmitter = SERVICE_NODE_MAPPING_ALONE;
        mappingCollector = SERVICE_NODE_MAPPING_ALONE;
        frequencyGovernor = cpufreq::GOVERNOR_USERSPACE;
        turboBoost = false;
        frequencyLowerBound = 0;
        frequencyUpperBound = 0;
        fastReconfiguration = false;
        migrateCollector = false;
        numSamples = 10;
        samplingInterval = 1000;
        underloadThresholdFarm = 80.0;
        overloadThresholdFarm = 90.0;
        underloadThresholdWorker = 80.0;
        overloadThresholdWorker = 90.0;
        requiredBandwidth = 0;
        maxBandwidthVariation = 5.0;
        requiredCompletionTime = 0;
        expectedTasksNumber = 0;
        powerBudget = 0;
        maxPredictionError = 10.0;
        voltageTableFile = "";
        observer = NULL;
    }
public:
    mammut::Mammut mammut; ///< The mammut modules handler.
    ContractType contractType; ///< The contract type that must be respected by the application
                               ///< [default = CONTRACT_NONE].
    StrategyMapping strategyMapping; ///< The mapping strategy [default = STRATEGY_MAPPING_LINEAR].
    StrategyHyperthreading strategyHyperthreading; ///< The hyperthreading strategy [default = STRATEGY_HT_NO].
    StrategyFrequencies strategyFrequencies; ///< The frequency strategy. It can be different from
                                             ///< STRATEGY_FREQUENCY_NO only if strategyMapping is
                                             ///< different from STRATEGY_MAPPING_NO [default = STRATEGY_FREQUENCY_NO].
    cpufreq::Governor frequencyGovernor; ///< The frequency governor (only used when
                                         ///< strategyFrequencies is STRATEGY_FREQUENCY_OS) [default = GOVERNOR_USERSPACE].
    bool turboBoost; ///< Flag to enable/disable cores turbo boosting [default = false].
    cpufreq::Frequency frequencyLowerBound; ///< The frequency lower bound (only if strategyFrequency is
                                            ///< STRATEGY_FREQUENCY_OS) [default = unused].
    cpufreq::Frequency frequencyUpperBound; ///< The frequency upper bound (only if strategyFrequency is
                                            ///< STRATEGY_FREQUENCY_OS) [default = unused].
    bool fastReconfiguration; ///< If true, before changing the number of workers the frequency will be set to
                              ///< maximum to reduce the latency of the reconfiguration. The frequency will be
                              ///< be set again to the correct value after the farm is restarted [default = false].
    StrategyUnusedVirtualCores strategyUnusedVirtualCores; ///< Strategy for virtual cores that are never used
                                                           ///< [default = STRATEGY_UNUSED_VC_NONE].
    StrategyUnusedVirtualCores strategyInactiveVirtualCores; ///< Strategy for virtual cores that become inactive
                                                             ///< after a workers reconfiguration
                                                             ///< [default = STRATEGY_UNUSED_VC_NONE].
    StrategyPrediction strategyPrediction; ///< Strategy to be used to predict power and performance values
                                           ///< [default = STRATEGY_PREDICTION_SIMPLE].
    ServiceNodeMapping mappingEmitter; ///< Emitter mapping [default = SERVICE_NODE_MAPPING_ALONE].
    ServiceNodeMapping mappingCollector; ///< Collector mapping [default = SERVICE_NODE_MAPPING_ALONE].
    bool migrateCollector; ///< If true, when a reconfiguration occur, the collector is migrated to a
                           ///< different virtual core (if needed) [default = false].
    uint32_t numSamples; ///< The number of samples used to take reconfiguration decisions [default = 10].
    uint32_t samplingInterval; ///< The length of the sampling interval (in milliseconds) for
                              ///< the data reading [default = 1000].
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
    double maxBandwidthVariation; ///< The allowed variation for bandwidth. The bandwidth will be kept
                                  ///< Between [B, B + x] where B is 'requiredBandwidth' and x
                                  ///< is the 'maxBandwidthVariation' percentage of B.
                                  ///< It is valid only if contractType is CONTRACT_BANDWIDTH [default = 5.0].
    uint requiredCompletionTime; ///< The required completion time for the application (in seconds). It is
                                 ///< valid only if contractType is CONTRACT_COMPLETION_TIME [default = unused].
    uint64_t expectedTasksNumber; ///< The number of task expected for this computation. It is
                                  ///< valid only if contractType is CONTRACT_COMPLETION_TIME [default = unused].
    double powerBudget;           ///< The maximum cores power to be used. It is
                                  ///< valid only if contractType is CONTRACT_POWER_BUDGET [default = unused].
    double maxPredictionError; ///< Maximum error percentage allowed for prediction [default = 10.0].
    std::string voltageTableFile; ///< The file containing the voltage table. It is mandatory when
                                  ///< strategyFrequencies is STRATEGY_FREQUENCY_YES or when contract is [default = unused].
    adp_ff_farm_observer* observer; ///< The observer object. It will be called every samplingInterval milliseconds
                                   ///< to monitor the adaptivity behaviour [default = NULL].

    /**
     * Creates the adaptivity parameters.
     * @param communicator The communicator used to instantiate the other modules.
     *        If NULL, the modules will be created as local modules.
     */
    AdaptivityParameters(Communicator* const communicator = NULL):
      mammut(communicator){
    	setDefault();
    }

    /**
     * Creates the adaptivity parameters.
     * @param xmlFileName The name containing the adaptivity parameters.
     * @param communicator The communicator used to instantiate the other modules.
     *        If NULL, the modules will be created as local modules.
     */
    AdaptivityParameters(const std::string& xmlFileName, Communicator* const communicator = NULL):
      mammut(communicator){
		setDefault();
		rapidxml::xml_document<> xmlContent;
		std::ifstream file(xmlFileName.c_str());
		if(!file.is_open()){
			throw std::runtime_error("Impossible to read xml file " + xmlFileName);
		}

		std::string fileContent;
		file.seekg(0, std::ios::end);
		fileContent.reserve(file.tellg());
		file.seekg(0, std::ios::beg);
		fileContent.assign((std::istreambuf_iterator<char>(file)),
							std::istreambuf_iterator<char>());
		char* fileContentChars = new char[fileContent.size() + 1];
		std::copy(fileContent.begin(), fileContent.end(), fileContentChars);
		fileContentChars[fileContent.size()] = '\0';
		xmlContent.parse<0>(fileContentChars);
		rapidxml::xml_node<> *root = xmlContent.first_node("adaptivityParameters");
		rapidxml::xml_node<> *node = NULL;

		node = root->first_node("contractType");
		if(node){
			contractType = (ContractType) utils::stringToInt(node->value());
		}

		node = root->first_node("strategyMapping");
		if(node){
			strategyMapping = (StrategyMapping) utils::stringToInt(node->value());
		}

        node = root->first_node("strategyHyperthreading");
        if(node){
            strategyHyperthreading = (StrategyHyperthreading) utils::stringToInt(node->value());
        }

		node = root->first_node("strategyFrequencies");
		if(node){
			strategyFrequencies = (StrategyFrequencies) utils::stringToInt(node->value());
		}

        node = root->first_node("strategyUnusedVirtualCores");
        if(node){
            strategyUnusedVirtualCores = (StrategyUnusedVirtualCores) utils::stringToInt(node->value());
        }

        node = root->first_node("strategyInactiveVirtualCores");
        if(node){
            strategyInactiveVirtualCores = (StrategyUnusedVirtualCores) utils::stringToInt(node->value());
        }

        node = root->first_node("strategyPrediction");
        if(node){
            strategyPrediction = (StrategyPrediction) utils::stringToInt(node->value());
        }

        node = root->first_node("mappingEmitter");
        if(node){
            mappingEmitter = (ServiceNodeMapping) utils::stringToInt(node->value());
        }

        node = root->first_node("mappingCollector");
        if(node){
            mappingCollector = (ServiceNodeMapping) utils::stringToInt(node->value());
        }

		node = root->first_node("frequencyGovernor");
		if(node){
			frequencyGovernor = (cpufreq::Governor) utils::stringToInt(node->value());
		}

		node = root->first_node("turboBoost");
		if(node){
			turboBoost = utils::stringToInt(node->value());
		}

		node = root->first_node("frequencyLowerBound");
		if(node){
			frequencyLowerBound = utils::stringToInt(node->value());
		}

		node = root->first_node("frequencyUpperBound");
		if(node){
			frequencyUpperBound = utils::stringToInt(node->value());
		}

		node = root->first_node("fastReconfiguration");
		if(node){
			fastReconfiguration = utils::stringToInt(node->value());
		}

		node = root->first_node("numSamples");
		if(node){
			numSamples = utils::stringToInt(node->value());
		}

		node = root->first_node("samplingInterval");
		if(node){
			samplingInterval = utils::stringToInt(node->value());
		}

		node = root->first_node("underloadThresholdFarm");
		if(node){
			underloadThresholdFarm = utils::stringToDouble(node->value());
		}

		node = root->first_node("overloadThresholdFarm");
		if(node){
			overloadThresholdFarm = utils::stringToDouble(node->value());
		}

		node = root->first_node("underloadThresholdWorker");
		if(node){
			underloadThresholdWorker = utils::stringToDouble(node->value());
		}

		node = root->first_node("overloadThresholdWorker");
		if(node){
			overloadThresholdWorker = utils::stringToDouble(node->value());
		}

		node = root->first_node("migrateCollector");
		if(node){
			migrateCollector = utils::stringToInt(node->value());
		}

		node = root->first_node("requiredBandwidth");
		if(node){
			requiredBandwidth = utils::stringToDouble(node->value());
		}

		node = root->first_node("maxBandwidthVariation");
		if(node){
			maxBandwidthVariation = utils::stringToDouble(node->value());
		}

		node = root->first_node("requiredCompletionTime");
        if(node){
            requiredCompletionTime = utils::stringToInt(node->value());
        }

        node = root->first_node("expectedTasksNumber");
        if(node){
            expectedTasksNumber = utils::stringToInt(node->value());
        }

        node = root->first_node("powerBudget");
        if(node){
            powerBudget = utils::stringToInt(node->value());
        }

        node = root->first_node("maxPredictionError");
        if(node){
            maxPredictionError = utils::stringToDouble(node->value());
        }

		node = root->first_node("voltageTableFile");
		if(node){
			voltageTableFile = node->value();
		}

		delete[] fileContentChars;
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
                if(requiredBandwidth < 0 || maxBandwidthVariation < 0 || maxBandwidthVariation > 100.0){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            case CONTRACT_PERF_COMPLETION_TIME:{
                if(!expectedTasksNumber || !requiredCompletionTime){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            case CONTRACT_POWER_BUDGET:{
                if(!powerBudget || (strategyFrequencies == STRATEGY_FREQUENCY_MIN_CORES) ||
                   strategyPrediction == STRATEGY_PREDICTION_SIMPLE){
                    return VALIDATION_WRONG_CONTRACT_PARAMETERS;
                }
            }break;
            default:{
                ;
            }break;
        }

        /** Validate voltage table. **/
        if(strategyFrequencies == STRATEGY_FREQUENCY_YES){
            if(voltageTableFile.empty() || !utils::existsFile(voltageTableFile)){
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

        /** Validate linear regression. **/
        if(strategyPrediction == STRATEGY_PREDICTION_REGRESSION_LINEAR){
            return VALIDATION_WRONG_REGRESSION_PARAMETERS;
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
