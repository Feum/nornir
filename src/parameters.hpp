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

/// Cores knob.
typedef enum{
    // The best solution is dinamically chosen.
    KNOB_CORES_AUTO = 0,

    // Doesn't changes the number of used cores.
    KNOB_CORES_NO,

    // Changes the number of threads used by the application.
    KNOB_CORES_CHANGE,

    // Instead of changing the number of threads it changes the mapping.
    KNOB_CORES_REMAP
}KnobCores;

/// Frequency knob.
typedef enum{
    // Chooses YES or NO according to the possibility of the system.
    KNOB_FREQUENCY_AUTO = 0,

    // Enables the possibility to dinamically change the frequencies.
    KNOB_FREQUENCY_YES,

    // Disables the possibility to dinamically change the frequencies.
    KNOB_FREQUENCY_NO,
}KnobFrequencies;

/// Mapping knob.
typedef enum{
    // Best mapping decisions will be at runtime.
    KNOB_MAPPING_AUTO = 0,

    // Mapping decisions will be performed by the operating system.
    KNOB_MAPPING_NO,

    // Tries to keep the threads as close as possible.
    KNOB_MAPPING_LINEAR,

    // Tries to make good use of the shared caches.
    // Particularly useful when threads have large working sets.
    KNOB_MAPPING_CACHE_EFFICIENT
}KnobMapping;

/// Hyperthreading knob.
typedef enum{
    // The runtime system will decide at runtime if use or not hyperthreading.
    KNOB_HT_AUTO = 0,

    // Hyperthreading is not used.
    KNOB_HT_NO = 0,

    // Hyperthreading is used since the beginning.
    KNOB_HT_YES_SOONER,

    // Hyperthreading is used only when we used all the physical cores.
    KNOB_HT_YES_LATER
}KnobHyperthreading;

/// Possible strategies to apply for unused virtual cores. For unused virtual
/// cores we mean those never used or those used only on some conditions.
typedef enum{
    // Automatically choose one of the other strategies.
    STRATEGY_UNUSED_VC_AUTO = 0,

    // Nothing is done on unused virtual cores.
    STRATEGY_UNUSED_VC_NONE,

    // Set the virtual cores to the lowest frequency (only
    // possible if all the other virtual cores on the same
    // domain are unused).
    STRATEGY_UNUSED_VC_LOWEST_FREQUENCY,

    // Turn off the virtual cores. They will not be anymore seen by the
    // operating system and it will not schedule anything on them.
    STRATEGY_UNUSED_VC_OFF
}StrategyUnusedVirtualCores;

/// Possible strategies to use to predict power and performance values.
typedef enum{
    // Applies multivariate linear regression.
    STRATEGY_PREDICTION_REGRESSION_LINEAR = 0,

    // Applies a simple analytical model.
    STRATEGY_PREDICTION_SIMPLE
}StrategyPrediction;

/// Possible ways to detect if the model has an high error.
typedef enum{
    // Constant prediction error, set by the user.
    STRATEGY_PREDICTION_ERROR_CONSTANT = 0,

    // Prediction error set equal to the maximum between
    // coefficient of variation and maximum specified error.
    STRATEGY_PREDICTION_ERROR_COEFFVAR
}StrategyPredictionError;

/// Service nodes (emitter or collector) mapping knob.
typedef enum{
    // The mapping is automatically chosen by the runtime system.
    KNOB_SNODE_MAPPING_AUTO = 0,

    // The service node is mapped on a physical core where no workers
    // are mapped.
    KNOB_SNODE_MAPPING_ALONE,

    // The service node is mapped on a physical core together with a worker.
    KNOB_SNODE_MAPPING_COLLAPSED,

    // The service node is mapped on an independent voltage domain and
    // is kept running at maximum performances (only available when
    // strategyFrequencies != STRATEGY_FREQUENCY_NO).
    KNOB_SNODE_MAPPING_PERFORMANCE
}KnobServiceNodeMapping;

/// Possible ways to smooth the values.
typedef enum{
    // Simple moving average
    STRATEGY_SMOOTHING_MOVING_AVERAGE = 0,

    // Exponential moving average
    STRATEGY_SMOOTHING_EXPONENTIAL
}StrategySmoothing;

/// Possible ways to change the smoothing factor at runtime.
typedef enum{
    // Constant smoothing factor
    STRATEGY_SMOOTHING_FACTOR_CONST = 0,

    // Changes the smoothing factor according to the variation
    // of the data.
    STRATEGY_SMOOTHING_FACTOR_DYNAMIC
}StrategySmoothingFactor;

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

// How to decide if we have to consider new configurations
typedef enum{
    // Number of samples.
    STRATEGY_PERSISTENCE_SAMPLES = 0,

    // Number of tasks.
    STRATEGY_PERSISTENCE_TASKS,

    // Coefficient of variation.
    STRATEGY_PERSISTENCE_VARIATION
}StrategyPersistence;

/// Possible parameters validation results.
typedef enum{
    // Parameters are ok.
    VALIDATION_OK = 0,

    // Frequency can be changed by the operating system and the flag
    // "constant_tsc" is not present on the CPU. Accordingly, since we
    // realy on getticks() to perform measurements, the amount of ticks
    // per second may change with the frequency, probably causing inaccuracy
    // in the measurements.
    VALIDATION_NO_CONSTANT_TSC,

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
}ParametersValidation;

/**
 * This class represents an XML tree.
 */
class XmlTree{
private:
    char* _fileContentChars;
    rapidxml::xml_node<>* _root;
public:
    /**
     * Reads the XML tree contained in 'fileName' and starting
     * with a root named 'rootName'.
     * @param fileName The name of the XML file.
     * @param rootName The name of the root.
     */
    XmlTree(const string& fileName, const string& rootName);

    ~XmlTree();

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a bool.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getBool(const char* valueName, bool& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an int.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getInt(const char* valueName, int& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an uint.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getUint(const char* valueName, uint& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an ulong.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getUlong(const char* valueName, ulong& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a double.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getDouble(const char* valueName, double& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a string.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getString(const char* valueName, string& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an enumeration.
     * @param T The type of the enumeration.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    template<typename T> void getEnum(const char* valueName, T& value);
};

#define SETVALUE(XML, TYPE, NAME) XML.get##TYPE(#NAME, NAME)

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

    void loadXml(const string& archFileName);
}ArchData;

/*!
 * \class AdaptivityParameters
 * \brief This class contains parameters for adaptivity choices.
 *
 * This class contains parameters for adaptivity choices.
 */
class Parameters{

private:
    friend class ManagerFarm;

    /**
     * Sets default parameters
     */
    void setDefault();

    /**
     * Sets the default values for parameters that depends
     * from others.
     */
    void setDefaultPost();

    /**
     * Checks if a governor is available on the current machine.
     * @param g The governor.
     * @return true if the governor is available, false otherwise.
     */
    bool isGovernorAvailable(Governor g);

    /**
     * Gets the frequencies available on this machine.
     * @return A vector containing the frequencies available on this machine.
     */
    vector<Frequency> getAvailableFrequencies();

    /**
     * Checks if it is possible to turn off the cores.
     * @return true if it is possible to turn off the cores,
     *         false otherwise.
     */
    bool isUnusedVcOffAvailable();

    /**
     * Checks if it is possible to set the frequency on this machine.
     * @return true if it is possible to set the frequency on this machine,
     *         false otherwise.
     */
    bool isFrequencySettable();

    /**
     * Checks if it is possible to set the lowest frequency on this machine.
     * @return true if it is possible to set the lowest frequency on this
     *         machine, false otherwise.
     */
    bool isLowestFrequencySettable();

    /**
     * Checks if it is possible to set the highest frequency on this machine.
     * @return true if it is possible to set the highest frequency on this
     *         machine, false otherwise.
     */
    bool isHighestFrequencySettable();

    /**
     * Checks if performance mapping has been required from one of the two
     * service nodes.
     * @return true if performance mapping has been required from one of
     *         the two service nodes, false otherwise.
     */
    bool serviceNodePerformance();

    /**
     * Validates the strategy required for unused virtual cores.
     * @param s The strategy required for unused virtual cores.
     * @return The result of the validation.
     */
    ParametersValidation validateUnusedVc(StrategyUnusedVirtualCores& s);

    /**
     * Validates the cores knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobCores();

    /**
     * Validates the frequency knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobFrequencies();

    /**
     * Validates the mapping knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobMapping();

    /**
     * Validates the service nodes knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobSnodeMapping();

    /**
     * Validates the hyperthreading knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobHt();

    /**
     * Validates the required contract.
     * @return The result of the validation.
     */
    ParametersValidation validateContract();

    /**
     * Loads the content of the parameters with the content
     * of an XML file.
     * @param fileName The name of the XML file.
     */
    void loadXml(const string& fileName);

public:
    // The mammut modules handler.
    Mammut mammut;

    // Architecture's specific data.
    ArchData archData;

    // The contract type that must be respected by the application
    // [default = CONTRACT_NONE].
    ContractType contractType;

    // The cores knob [default = KNOB_CORES_AUTO].
    KnobCores knobCores;

    // The frequency knob. It can be KNOB_FREQUENCY_YES
    // only if knobMapping is different from KNOB_MAPPING_NO
    // [default = KNOB_FREQUENCY_AUTO].
    KnobFrequencies knobFrequencies;

    //  The mapping knob [default = KNOB_MAPPING_AUTO].
    KnobMapping knobMapping;

    // Emitter mapping knob [default = KNOB_SNODE_MAPPING_AUTO].
    KnobServiceNodeMapping knobMappingEmitter;

    // Collector mapping knob [default = KNOB_SNODE_MAPPING_AUTO].
    KnobServiceNodeMapping knobMappingCollector;

    // The hyperthreading knob [default = KNOB_HT_AUTO].
    KnobHyperthreading knobHyperthreading;

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
    Parameters(Communicator* const communicator = NULL);

    /**
     * Creates the adaptivity parameters.
     * @param paramFileName The name of the XML file containing the adaptivity
     *        parameters.
     * @param archFileName The name of the XML file containing the
     *        architectural data.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    Parameters(const string& paramFileName,
               const string& archFileName,
               Communicator* const communicator = NULL);

    /**
     * Destroyes these parameters.
     */
    ~Parameters();

    /**
     * Validates these parameters.
     * @return The validation result.
     */
    ParametersValidation validate();
};

}

#endif /* ADAPTIVE_FASTFLOW_PARAMETERS_HPP_ */
