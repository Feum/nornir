/*
 * parameters.hpp
 *
 * Created on: 23/03/2015
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
 *
 * Strategies used/developed for different papers:
 * - PDP2015:
 *         strategySelection = STRATEGY_SELECTION_ANALYTICAL
 *         knobMappingEnabled = false
 *         (Actually was a simpler selection strategy than the current one).
 *     - PPL2016:
 *         strategySelection = STRATEGY_SELECTION_ANALYTICAL
 *         knobMappingEnabled = false
 * - PDP2016:
 *         strategySelection =  STRATEGY_SELECTION_LEARNING
 *         strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_AMDAHL
 *         strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR
 *         strategyExploration = (Many of them)
 *         knobMappingEnabled = false
 * - TACO2016:
 *         strategySelection = STRATEGY_SELECTION_LEARNING
 *         strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_USL
 *         strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR
 *         strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY
 *         strategyExploration = STRATEGY_EXPLORATION_HALTON
 *         knobMappingEnabled = true
 *
 *         In addition to that, we also developed and compared with:
 *             strategySelection = STRATEGY_SELECTION_MISHRA
 *             strategySelection = STRATEGY_SELECTION_FULLSEARCH
 *             strategySelection = STRATEGY_SELECTION_LIMARTINEZ
 *         and also with
 *             knobMappingEnabled = false
 *
 */


#ifndef NORNIR_PARAMETERS_HPP_
#define NORNIR_PARAMETERS_HPP_

#include "external/rapidXml/rapidxml.hpp"

#include <cmath>
#include <fstream>
#include "external/mammut/mammut/utils.hpp"
#include "external/mammut/mammut/mammut.hpp"

namespace nornir{

#define MANAGER_VIRTUAL_CORE (VirtualCoreId) 0

class Observer;

// Possible knobs
typedef enum{
    KNOB_TYPE_VIRTUAL_CORES = 0, // Number of contexts to be used.
    KNOB_TYPE_HYPERTHREADING, // Number of contexts to be used on each physical core.
    KNOB_TYPE_MAPPING, // Mapping of threads on physical cores.
    KNOB_TYPE_FREQUENCY,
    KNOB_TYPE_NUM  // <---- This must always be the last value
}KnobType;

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

    // The maximum bandwidth is requested. No considerations about power consumption
    // will be done.
    CONTRACT_PERF_MAX,

    // A specific maximum cores power is requested. Under this constraint, the
    // configuration with best performance is chosen.
    CONTRACT_POWER_BUDGET
}ContractType;

/// Communication queues blocking/nonblocking.
typedef enum{
    // Non blocking queue.
    TRIGGER_Q_BLOCKING_NO = 0,

    // Blocking queue.
    TRIGGER_Q_BLOCKING_YES,

    // Automatically choose between blocking and nonblocking.
    TRIGGER_Q_BLOCKING_AUTO
}TriggerConfQBlocking;

/// Possible strategies to apply for unused virtual cores.
typedef enum{
    // Automatically choose one of the other strategies.
    STRATEGY_UNUSED_VC_AUTO = 0,

    // The unused virtual cores will run at the same
    // frequency of the used ones.
    STRATEGY_UNUSED_VC_SAME,

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

// Possible strategies for the selection of the best configuration.
typedef enum{
    // Applies an online learning algorithm
    STRATEGY_SELECTION_LEARNING = 0,

    // Applies a simple analytical model.
    STRATEGY_SELECTION_ANALYTICAL,

    // Tries all the configurations in order to find the best one.
    STRATEGY_SELECTION_FULLSEARCH,

    // Applies the algorithm described in:
    // "Dynamic Power-Performance Adaptation of Parallel Computation
    // on Chip Multiprocessors" - Jian Li and Jose F. Martınez
    STRATEGY_SELECTION_LIMARTINEZ,

    // Applies the algorithm described in:
    // "A Probabilistic Graphical Model-based Approach for Minimizing
    // Energy Under Performance Constraints" - Mishra, Nikita and Zhang, Huazhe
    // and Lafferty, John D. and Hoffmann, Henry
    STRATEGY_SELECTION_MISHRA,

    STRATEGY_SELECTION_NUM // <- Must always be the last.
}StrategySelection;

// Possible prediction strategies for performance. Can only be specified if the
// selection strategy is "LEARNING".
typedef enum{
    STRATEGY_PREDICTION_PERFORMANCE_AMDAHL = 0,
    STRATEGY_PREDICTION_PERFORMANCE_USL,
    STRATEGY_PREDICTION_PERFORMANCE_USLP, // <- More precise than USL but needs one additional calibration point.
    STRATEGY_PREDICTION_PERFORMANCE_MISHRA,
    STRATEGY_PREDICTION_PERFORMANCE_NUM // <- This must always be the last.
}StrategyPredictionPerformance;

// Possible prediction strategies for power consumption. Can only be specified
// if the selection strategy is "LEARNING".
typedef enum{
    STRATEGY_PREDICTION_POWER_LINEAR = 0,
    STRATEGY_PREDICTION_POWER_MISHRA,
    STRATEGY_PREDICTION_POWER_NUM // <- This must always be the last.
}StrategyPredictionPower;

/// Possible ways to select the calibration points. Can only be specified if
/// the selection strategy is "LEARNING".
typedef enum{
    // Random choice of the points.
    STRATEGY_EXPLORATION_RANDOM = 0,
    // Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992)
    STRATEGY_EXPLORATION_NIEDERREITER,
    // Antonov, Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980)
    STRATEGY_EXPLORATION_SOBOL,
    // J.H. Halton, Numerische Mathematik 2, 84-90 (1960) and B. Vandewoestyne
    // R. Cools Computational and Applied Mathematics 189, 1&2, 341-361 (2006)
    STRATEGY_EXPLORATION_HALTON,
    STRATEGY_EXPLORATION_HALTON_REVERSE,
}StrategyExploration;

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

    // Coefficient of variation.
    STRATEGY_PERSISTENCE_VARIATION
}StrategyPersistence;

/// Possible parameters validation results.
typedef enum{
    // Parameters are ok.
    VALIDATION_OK = 0,

    // Generic error
    VALIDATION_NO,

    // Some of the knobs enabled by the user is not supported
    // by the selection algorithm he chosen.
    VALIDATION_UNSUPPORTED_KNOBS,

    // On this architecture is not possible to manually change
    // the frequency.
    VALIDATION_NO_MANUAL_DVFS,

    // On this architecture fast reconfiguration is not available.
    VALIDATION_NO_FAST_RECONF,

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

    // Blocking threshold needs to be specified.
    VALIDATION_NO_BLOCKING_THRESHOLD,

    // Parameters for Mishra predictors not specified.
    VALIDATION_NO_MISHRA_PARAMETERS,
}ParametersValidation;

/**
 * This class represents an XML tree.
 */
class XmlTree{
private:
    char* _fileContentChars;
    rapidxml::xml_node<>* _root;

    void init(const std::string& content, const std::string& rootName);
public:
    /**
     * Reads the XML tree contained in 'fileName' and starting
     * with a root named 'rootName'.
     * @param fileName The name of the XML file.
     * @param rootName The name of the root.
     */
    XmlTree(const std::string& fileName, const std::string& rootName);

    ~XmlTree();


    /**
     * Returns a node with name valueName.
     * @param valueName The name of the node. It can be of the form
     *        field.subfield.subsubfield. ... (similar to a C struct).
     * @return The node with name valueName.
     */
    rapidxml::xml_node<>* getNode(const char* valueName) const;

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
    void getString(const char* valueName, std::string& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an array of uint
     * (separated by ':').
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getArrayUint(const char* valueName, std::vector<uint>& value);

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

    mammut::cpufreq::VoltageTable voltageTable;

    ArchData():ticksPerNs(0),
               monitoringCost(0){;}

    void loadXml(const std::string& archFileName);
}ArchData;

typedef struct{
    /**
     * The name of this application (it must matches with one of the names
     * in applicationNames file).
     */
    std::string applicationName;
    /**
     * File containing the names of the applications, in the same
     * order of the columns in bandwidth and power data.
     */
    std::string namesData;
    /**
     * File containing the data about the bandwidth. One column per application,
     * one row per configuration. Data doesn't need to be normalized.
     */
    std::string bandwidthData;
    /**
     * File containing the data about the power. One column per application,
     * one row per configuration. Data doesn't need to be normalized.
     */
    std::string powerData;
    /**
     * Number of samples to be used [default = 20].
     */
    uint numSamples;
}MishraParameters;

typedef struct{
    /**
     * If true, the interpreter will ensure that instructions belonging
     * to different stream elements will be processed in the same order
     * they are received [default = false].
     */
    bool orderedProcessing;

    /**
     * If true, the interpreter will ensure that the output will be produced
     * in the same order of the received data [default = false]. If
     * orderedProcessing is true, this property is automatically guaranteed.
     */
    bool orderedOutput;

    /**
     * Maximum number of graphs to keep in the system [default = 1000].
     */
    uint maxGraphs;

    /**
     * Maximum number of interpreters to be used [default = #Physical cores
     * in the system - 2].
     */
    uint maxInterpreters;
}DataflowParameters;

/*!
 * \class AdaptivityParameters
 * \brief This class contains parameters for adaptivity choices.
 *
 * This class contains parameters for adaptivity choices.
 */
class Parameters{
private:
    bool _knobEnabled[KNOB_TYPE_NUM];

    /**
     * Sets default parameters
     */
    void setDefault();

    /**
     * Computes the sampling interval such to have a low
     * performance overhead.
     * @return The sampling interval.
     */
    uint getLowOverheadSamplingInterval() const;

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
    bool isGovernorAvailable(mammut::cpufreq::Governor g);

    /**
     * Gets the frequencies available on this machine.
     * @return A vector containing the frequencies available on this machine.
     */
    std::vector<mammut::cpufreq::Frequency> getAvailableFrequencies();

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
     * Validates the strategy required for unused virtual cores.
     * @param s The strategy required for unused virtual cores.
     * @return The result of the validation.
     */
    ParametersValidation validateUnusedVc(StrategyUnusedVirtualCores& s);

    /**
     * Validates the frequency knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobFrequencies();

    /**
     * Validates the triggers.
     * @return The result of the validation.
     */
    ParametersValidation validateTriggers();

    /**
     * Validates the required contract.
     * @return The result of the validation.
     */
    ParametersValidation validateContract();

    /**                                                                                                                                                
     * Validates the selector.
     * @return The result of the validation.
     */
    ParametersValidation validateSelector();

    /**
     * Loads the content of the parameters with the content
     * of an XML file.
     * @param fileName The name of the XML file.
     */
    void loadXml(const std::string& fileName);

public:
    // The mammut modules handler.
    mammut::Mammut mammut;

    // Architecture's specific data.
    ArchData archData;

    // The contract type that must be respected by the application
    // [default = CONTRACT_NONE].
    ContractType contractType;

    // The Q blocking knob [default = KNOB_Q_BLOCKING_NO].
    TriggerConfQBlocking triggerQBlocking;

    // Strategy for unused virtual cores [default = STRATEGY_UNUSED_VC_SAME].
    StrategyUnusedVirtualCores strategyUnusedVirtualCores;

    // Strategy to be used to select the best configuration according
    // to the requirements [default = STRATEGY_SELECTION_LEARNING].
    StrategySelection strategySelection;

    // Strategy to be used to predict performance values.
    // Only valid when strategySelection is LEARNING
    // [default = STRATEGY_PREDICTION_PERFORMANCE_USL].
    StrategyPredictionPerformance strategyPredictionPerformance;

    // Strategy to be used to predict power values.
    // Only valid when strategySelection is LEARNING
    // [default = STRATEGY_PREDICTION_POWER_LINEAR].
    StrategyPredictionPower strategyPredictionPower;

    // Strategy to be used to select the points to be explored during
    // calibration. Only valid when strategySelection is LEARNING
    // [default = STRATEGY_EXPLORATION_HALTON].
    StrategyExploration strategyExploration;

    // Smoothing strategy [default = STRATEGY_SMOOTHING_EXPONENTIAL].
    StrategySmoothing strategySmoothing;

    // Polling strategy [default = STRATEGY_POLLING_SLEEP_SMALL].
    StrategyPolling strategyPolling;

    // Persistence strategy [default = STRATEGY_PERSISTENCE_SAMPLES].
    StrategyPersistence strategyPersistence;

    // Flag to enable/disable cores knobs [default = true].
    bool knobCoresEnabled;

    // Flag to enable/disable mapping knob [default = true].
    bool knobMappingEnabled;

    // Flag to enable/disable frequency knob [default = true].
    bool knobFrequencyEnabled;

    // Flag to enable/disable hyperthreading knob [default = false].
    bool knobHyperthreadingEnabled;

    // Parameters for Mishra predictor.
    MishraParameters mishra;

    // Flag to enable/disable cores turbo boosting [default = false].
    bool turboBoost;

    // If true, before changing the number of workers the frequency will be
    // set to maximum to reduce the latency of the reconfiguration. The
    // frequency will be be set again to the correct value after the farm
    // is restarted [default = true].
    bool fastReconfiguration;

    // If true, when a reconfiguration occur, the collector is migrated to a
    // different virtual core (if needed) [default = false].
    bool migrateCollector;

    // The smoothing factor. It's meaning changes according to the smoothing
    // strategy adopted. If 0, default value will be set.
    // [default = 10 for moving average, 0.1 for exponential].
    double smoothingFactor;

    // The persistence value. It's meaning changes according to the persistence
    // strategy adopted. If 0, default value will be set.
    // [default = 1 for samples, 5 for variation].
    double persistenceValue;

    // The length of the sampling interval (in milliseconds) for the data
    // reading during calibration phase. If 0, it will be automatically computed
    // such to have a low performance overhead [default = 100].
    uint32_t samplingIntervalCalibration;

    // The length of the sampling interval (in milliseconds) for the data
    // reading during steady phase. If 0, it will be automatically computed
    // such to have a low performance overhead [default = 1000].
    uint32_t samplingIntervalSteady;

    // If the application stays in the steady phase for at least
    // steadyThreshold samples, then we consider the application to be stedy
    // [default = 4].
    uint32_t steadyThreshold;

    /// The minimum number of tasks in a worker sample. If 0, no minimum.
    /// [default = 0].
    uint minTasksPerSample;

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

    /// If true, this is an application when the workers works synchronously,
    /// i.e. the emitter works as a barrier between two successive iterations
    /// of the workers. In this case, the emitter always broadcasts tasks
    /// to all the workers [default = false].
    bool synchronousWorkers;

    // The maximum cores power to be used. It is
    // valid only if contractType is CONTRACT_POWER_BUDGET [default = unused].
    double powerBudget;

    // Maximum calibration time (milliseconds). 0 is no limit.
    // We will keep calibrating until the error is higher than the
    // max*PredictionError AND calibration time is lower than the
    // maxCalibrationTime  AND the number of visited configurations
    // is lower than maxCalibrationConfigurations [default = 0.0].
    double maxCalibrationTime;

    // Maximum number of configurations to explore during calibration
    // phase. 0 is no limit. We will keep calibrating until the error is
    // higher than the max*PredictionError AND calibration time is lower
    // than the maxCalibrationTime AND the number of visited configurations
    // is lower than maxCalibrationConfigurations [default = 0.0].
    uint maxCalibrationConfigurations;

    // Maximum error percentage allowed for performance prediction.
    // [default = 10.0].
    double maxPerformancePredictionError;

    // Maximum error percentage allowed for power consumption prediction.
    // [default = 5.0].
    double maxPowerPredictionError;

    // The aging value for the linear regression predictor. If it has a value
    // of n, then we will only consider the last n configurations we collected
    // in order to perform predictions. If 0, no aging will be applied and all
    // the previous samples will be considered [default = 0].
    uint regressionAging;

    // The maximum percentage of monitoring overhead, in the range (0, 100).
    // [default = 1.0].
    double maxMonitoringOverhead;

    // Idle threshold (in microseconds) to switch from non blocking to blocking
    // runtime support. If the current idle time goes above this value,
    // and if the runtime has been configured to do so, it will switch
    // from non blocking to blocking support (and viceversa). If -1.0,
    // the runtime will never switch [default = -1.0].
    double thresholdQBlocking;

    // Max number of samples that can be violated (either accuracy or contract
    // violations) [default = 0].
    uint tolerableSamples;

    // Maximum size for the internal queues of the farm.
    // 0 corresponds to infinite size [default = 1].
    ulong qSize;

    // If different from zero:
    //      - If we have a PERF_* contract, we will search for a solution which
    //        is the conservativeValue% more performing than the requirement.
    //      - If we have a POWER_* contract, we will search for a solution
    //        which consumes conservativeValue% less power than the requirement.
    // This is done to amortize fluctuations. [default = 0.0]
    double conservativeValue;

    // A vector containing the number of cores not allowed to be used.
    // E.g. if it contains the number 3, then the runtime will
    // never use  3 cores. It can only be specified when
    // knobCoresEnabled is true [default = empty].
    std::vector<uint> disallowedNumCores;

    // If true, the manager will run on a physical core by itself.
    // If false, it will run on the same physical core used by the emitter
    // of the farm [default = false].
    bool isolateManager;

    // If true, computes the statistics about the cost
    // of the reconfigurations [default = false].
    bool statsReconfiguration;

    // Parameters for dataflow applications.
    DataflowParameters dataflow;

    // The observer object. It will be called every samplingInterval
    // milliseconds to monitor the adaptivity behaviour [default = NULL].
    Observer* observer;

    /**
     * Creates the adaptivity parameters.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    Parameters(mammut::Communicator* const communicator = NULL);

    /**
     * Creates the adaptivity parameters.
     * @param paramFileName The name of the XML file containing the adaptivity
     *        parameters.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    Parameters(const std::string& paramFileName,
               mammut::Communicator* const communicator = NULL);

    /**
     * Destroyes these parameters.
     */
    ~Parameters();

    /**
     * Validates these parameters.
     * @return The validation result.
     */
    ParametersValidation validate();

    /**
     * Check if a specific knob is enabled.
     * ATTENTION: Can only be called after validate().
     * @param k The knob to check.
     * @return true if the specified knob is enabled, false otherwise.
     */
    bool isKnobEnabled(KnobType k) const;
};

}

#endif /* NORNIR_PARAMETERS_HPP_ */
