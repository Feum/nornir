/*
 * manager.hpp
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

/*! \file manager.hpp
 * \brief Implementation of an adaptive fastflow farm.
 *
 * To let an existing fastflow farm-based adaptive, follow these steps:
 *  1. Emitter, Workers and Collector of the farm must extend
 *     adpff::AdaptiveNode instead of ff::ff_node
 *  2. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyWorkersChange virtual
 *     method.
 */

#ifndef ADAPTIVE_FASTFLOW_FARM_HPP_
#define ADAPTIVE_FASTFLOW_FARM_HPP_

#include "parameters.hpp"
#include "predictors.hpp"
#include "node.hpp"
#include "utils.hpp"

#include <ff/farm.hpp>
#include <mammut/module.hpp>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

#include <cmath>
#include <iostream>
#include <limits>

namespace adpff{

class Parameters;
class ManagerFarm;

using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

struct MonitoredSample;

/*!
 * \internal
 * \struct FarmConfiguration
 * \brief Represents a possible farm configuration.
 *
 * This struct represent a possible farm configuration.
 */
typedef struct FarmConfiguration{
    uint numWorkers;
    Frequency frequency;

    FarmConfiguration():numWorkers(1), frequency(1){;}
    FarmConfiguration(uint numWorkers, Frequency frequency = 1):
                     numWorkers(numWorkers), frequency(frequency){;}
}FarmConfiguration;

/*!
 * \class AdaptivityManagerFarm
 * \brief This class manages the adaptivity in farm based computations.
 *
 * This class manages the adaptivity in farm based computations.
 */
class ManagerFarm: public Thread{
    friend class PredictorSimple;
    friend class PredictorLinearRegression;
    friend class RegressionData;
    friend class RegressionDataServiceTime;
    friend class RegressionDataPower;
    friend class Calibrator;
    friend class CalibratorLowDiscrepancy;

public:
    /**
     * Creates a farm adaptivity manager.
     * @param farm The farm to be managed.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFarm(ff_farm<>* farm, Parameters adaptivityParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFarm();

    /**
     * Function executed by this thread.
     */
    void run();
private:
    // The managed farm.
    ff_farm<>* _farm;

    // The parameters used to take management decisions.
    Parameters _p;

    // Starting time of the manager.
    uint _startTimeMs;

    // The cpufreq module.
    CpuFreq* _cpufreq;

    // The energy module.
    Energy* _energy;

    // The task module.
    TasksManager* _task;

    // The topology module.
    Topology* _topology;

    // Number of CPUs
    uint _numCpus;

    // Number of physical cores.
    uint _numPhysicalCores;

    // Number of physical cores per CPU.
    uint _numPhysicalCoresPerCpu;

    // Number of virtual cores per physical core.
    uint _numVirtualCoresPerPhysicalCore;

    // The emitter (if present).
    AdaptiveNode* _emitter;

    // The collector (if present).
    AdaptiveNode* _collector;

    // The currently running workers.
    vector<AdaptiveNode*> _activeWorkers;

    // Workers that can run but are not currently running.
    vector<AdaptiveNode*> _inactiveWorkers;

    // The maximum number of workers that can be activated by the manager.
    size_t _maxNumWorkers;

    // If true, the user requested sensitivity for emitter and the
    // request has been satisfied.
    bool _emitterSensitivitySatisfied;

    // If true, the user requested sensitivity for collector and the
    // request has been satisfied.
    bool _collectorSensitivitySatisfied;

    // The current configuration of the farm.
    FarmConfiguration _currentConfiguration;

    // CPUs currently used by farm nodes.
    vector<CpuId> _usedCpus;

    // CPUs not used by farm nodes.
    vector<CpuId> _unusedCpus;

    // The available virtual cores, sorted according to  the mapping strategy.
    const vector<VirtualCore*> _availableVirtualCores;

    // The virtual cores where the active workers are running.
    vector<VirtualCore*> _activeWorkersVirtualCores;

    ///< The virtual cores where the inactive workers are running.
    vector<VirtualCore*> _inactiveWorkersVirtualCores;

    // Virtual cores not used by the farm nodes.
    vector<VirtualCore*> _unusedVirtualCores;

    // The virtual core where the emitter (if present) is running.
    VirtualCore* _emitterVirtualCore;

    // The virtual core where the collector (if present) is running.
    VirtualCore* _collectorVirtualCore;

    // The domains on which frequency scaling is applied.
    vector<Domain*> _scalableDomains;

    // The voltage table.
    VoltageTable _voltageTable;

    // The available frequencies on this machine.
    vector<Frequency> _availableFrequencies;

    // Monitored samples;
    Smoother<MonitoredSample>* _samples;

    // The number of tasks processed since the last reconfiguration.
    double _totalTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the number of tasks
    // that still needs to be processed by the application.
    uint64_t _remainingTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the deadline of
    // the application.
    time_t _deadline;

    // Milliseconds timestamp of the last store of a sample.
    double _lastStoredSampleMs;

    // The calibrator of the predictors.
    Calibrator* _calibrator;

    // The predictor of the primary value.
    Predictor* _primaryPredictor;

    // The predictor of the secondary value.
    Predictor* _secondaryPredictor;

    // The prediction done for the primary value for the chosen configuration.
    double _primaryPrediction;

    // The prediction done for the secondary value for the chosen configuration.
    double _secondaryPrediction;

#ifdef DEBUG_FARM
    ofstream samplesFile;
#endif

    /**
     * Returns true if the number of workers must be reconfigured.
     * @return true if the number of workers must be reconfigured,
     *         false otherwise
     */
    bool reconfigureWorkers() const;

    /**
     * Returns true if the frequencies must be reconfigured.
     * @return true if the frequencies must be reconfigured,
     *         false otherwise
     */
    bool reconfigureFrequency() const;

    /**
     * Returns the number of dimensions of a configuration,
     * i.e. the number of different decisions that can
     * be taken at runtime.
     */
    uint getConfigurationDimension() const;

    /**
     * If possible, finds a set of physical cores belonging to domains
     * different from those of virtual cores in 'virtualCores' vector.
     * @param virtualCores A vector of virtual cores.
     * @return A set of physical cores that can always run at the highest
     *         frequency.
     */
    vector<PhysicalCore*> getSepDomainsPhyCores(const vector<VirtualCore*>&
                                                      virtualCores) const;

    /**
     * Set a specified domain to the highest frequency.
     * @param domain The domain.
     */
    void setDomainToHighestFrequency(const Domain* domain);

    /**
     * Computes the available virtual cores, sorting them according to
     * the specified mapping strategy.
     * @return The available virtual cores, sorted according to the
     *         specified mapping strategy.
     */
    vector<VirtualCore*> getAvailableVirtualCores();

    /**
     * Returns the number of scalable service nodes.
     * @return The number of scalable service nodes.
     */
    uint numScalableServiceNodes();

    /**
     * Manages mapping of emitter and collector.
     */
    void manageServiceNodesPerformance();

    /**
     * Generates mapping indexes. They are indexes to be used on
     * _availableVirtualCores vector to get the corresponding virtual core
     * where a specific node must be mapped.
     * @param emitterIndex The index of the emitter.
     * @param firstWorkerIndex The index of the first worker
     *                         (the others follow).
     * @param collectorIndex The index of the collector (if present).
     */
    void getMappingIndexes(size_t& emitterIndex,
                           size_t& firstWorkerIndex,
                           size_t& collectorIndex);

    /**
     * Computes the virtual cores where the nodes must be mapped
     * and pins these nodes on the virtual cores.
     */
    void mapNodesToVirtualCores();

    /**
     * Applies the OFF strategy for a specified set of virtual cores.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategyOff(const vector<VirtualCore*>& unusedVc);

    /**
     * Applies the LOWEST_FREQUENCY strategy for a specified set of
     * virtual cores.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategyLowestFreq(const vector<VirtualCore*>& vc);

    /**
     * Applies a specific strategy for a specified set of virtual cores.
     * @param strategyUnused The unused strategy.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategy(StrategyUnusedVirtualCores strategyUnused,
                               const vector<VirtualCore*>& vc);
    /**
     * Apply the strategies for inactive and unused virtual cores.
     */
    void applyUnusedVCStrategy();

    /**
     * Updates the scalable domains vector.
     */
    void updateScalableDomains();

    /**
     * Set a specific P-state for the virtual cores used by
     * the current active workers, emitter (if not sensitive) and
     * collector (if not sensitive).
     * @param frequency The frequency to be set.
     */
    void updatePstate(Frequency frequency);

    /**
     * Map the nodes to virtual cores and
     * prepares frequencies and governors for running.
     */
    void mapAndSetFrequencies();

    /**
     * Returns the maximum primary prediction error.
     * @return The maximum primary prediction error.
     */
    double getMaxPredictionErrorPrimary() const;

    /**
     * Returns the maximum secondary prediction error.
     * @return The maximum secondary prediction error.
     */
    double getMaxPredictionErrorSecondary() const;

    /**
     * Returns the primary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The primary value of a sample according to
     * the required contract.
     */
    double getPrimaryValue(const MonitoredSample& sample) const;

    /**
     * Returns the secondary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The secondary value of a sample according to
     * the required contract.
     */
    double getSecondaryValue(const MonitoredSample& sample) const;
    /**
     * Returns the primary value according to the required contract.
     * @return The primary value according to the required contract.
     */
    double getPrimaryValue() const;

    /**
     * Returns the secondary value according to the required contract.
     * @return The secondary value according to the required contract.
     */
    double getSecondaryValue() const;

    /**
     * Checks if the contract is violated.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated() const;

    /**
     * Checks if a specific primary value respects the required contract.
     * @param value The value to be checked.
     * @param tolerance The percentage of tolerance allowed for the check
     */
    bool isFeasiblePrimaryValue(double value, double tolerance = 0) const;

    /**
     * Returns the voltage at a specific configuration.
     * @param configuration The configuration.
     * @return The voltage at that configuration.
     */
    double getVoltage(const FarmConfiguration& configuration) const;

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

    /**
     * Computes the new configuration of the farm after a contract violation.
     * @return The new configuration.
     */
    FarmConfiguration getNewConfiguration();
    /**
     * Updates the currently used and unused CPUs.
     */
    void updateUsedCpus();

    /**
     * Changes the active and inactive nodes according to the new configuration.
     * @param configuration The new configuration.
     */
    void changeActiveNodes(FarmConfiguration& configuration);

    /**
     * Prepares the nodes to freeze.
     */
    void prepareToFreeze();

    /**
     * Freezes all the nodes.
     */
    void freeze();

    /**
     * Notifies a change in the number of workers to all the nodes.
     * @param numWorkers The new number of workers.
     */
    void notifyNewConfiguration(uint numWorkers);

    /**
     * Prepares the nodes to run.
     * @param numWorkers The new number of workers.
     */
    void prepareToRun(uint numWorkers);

    /**
     * Runs the farm with a new number of workers.
     * @param numWorkers The new number of workers.
     */
    void run(uint numWorkers);

    /**
     * Checks if the application terminated.
     */
    bool terminated();

    /**
     * Changes the current farm configuration.
     * @param configuration The new configuration.
     */
    void changeConfiguration(FarmConfiguration configuration);

    /**
     * Send data to observer.
     **/
    void observe();

    /**
     * Asks the workers for their samples.
     */
    void askForWorkersSamples();

    /**
     * Obtain workers samples.
     * @param sample A worker sample. It will be filled by this call with the
     *               global data of the farm.
     */
    void getWorkersSamples(WorkerSample& sample);

    /**
     * Store a new sample.
     **/
    void storeNewSample();

    /**
     * Returns true if the manager doesn't have still to check for a new
     * configuration.
     * @return True if the manager doesn't have still to check for a new
     * configuration.
     */
    bool persist() const;

    /**
     * Initializes the calibrator.
     */
    void initCalibrator();

    /**
     * Initializes the predictors
     */
    void initPredictors();
};

/*!
 * This class can be used to obtain statistics about reconfigurations
 * performed by the manager.
 * It can be extended by a user defined class to customize action to take
 * for each observed statistic.
 */
class Observer{
    friend class ManagerFarm;
private:
    ofstream _statsFile;
    ofstream _calibrationFile;
    ofstream _summaryFile;
    unsigned int _startMonitoringMs;
    double _totalWatts;
    double _totalBw;
    unsigned long _numSamples;

    double calibrationDurationToPerc(const CalibrationStats& cs,
                                     uint durationMs);
public:
    Observer(string statsFile = "stats.csv",
             string calibrationFile = "calibration.csv",
             string summaryFile = "summary.csv");

    virtual ~Observer();

    virtual void observe(unsigned int timeStamp,
                         size_t workers,
                         Frequency frequency,
                         const VirtualCore* emitterVirtualCore,
                         const vector<VirtualCore*>& workersVirtualCore,
                         const VirtualCore* collectorVirtualCore,
                         double currentBandwidth,
                         double smoothedBandwidth,
                         double coeffVarBandwidth,
                         double smoothedUtilization,
                         JoulesCpu smoothedWatts);

    virtual void calibrationStats(const vector<CalibrationStats>&
                                  calibrationStats,
                                  uint durationMs);

    virtual void summaryStats(const vector<CalibrationStats>&
                              calibrationStats,
                              uint durationMs);
};

inline ostream& operator<<(ostream& os, const FarmConfiguration& obj){
    os << "[";
    os << obj.numWorkers << ", ";
    os << obj.frequency;
    os << "]";
    return os;
}

inline bool operator==(const FarmConfiguration& lhs,
                       const FarmConfiguration& rhs){
    return lhs.numWorkers == rhs.numWorkers &&
           lhs.frequency == rhs.frequency;
}

inline bool operator!=(const FarmConfiguration& lhs,
                       const FarmConfiguration& rhs){
    return !operator==(lhs,rhs);
}

inline bool operator<(const FarmConfiguration& lhs,
                      const FarmConfiguration& rhs){
    if(lhs.numWorkers < rhs.numWorkers){
        return true;
    }
    if(lhs.numWorkers > rhs.numWorkers){
        return false;
    }

    if(lhs.frequency < rhs.frequency){
        return true;
    }
    return false;
}

inline bool operator>(const FarmConfiguration& lhs,
                      const FarmConfiguration& rhs){
    return operator< (rhs,lhs);
}

inline bool operator<=(const FarmConfiguration& lhs,
                       const FarmConfiguration& rhs){
    return !operator> (lhs,rhs);
}

inline bool operator>=(const FarmConfiguration& lhs,
                       const FarmConfiguration& rhs){
    return !operator< (lhs,rhs);
}

typedef struct MonitoredSample{
    JoulesCpu watts; ///< Watts consumed by all the CPUs.
    double bandwidth; ///< Bandwidth of the entire farm.
    double utilization; ///< Utilization of the entire farm.
    double latency; ///< Average latency of a worker (nanoseconds).

    MonitoredSample():bandwidth(0), utilization(0), latency(0){;}

    void swap(MonitoredSample& x){
        using std::swap;

        swap(watts, x.watts);
        swap(bandwidth, x.bandwidth);
        swap(utilization, x.utilization);
        swap(latency, x.latency);
    }

    MonitoredSample& operator=(MonitoredSample rhs){
        swap(rhs);
        return *this;
    }

    MonitoredSample& operator+=(const MonitoredSample& rhs){
        watts += rhs.watts;
        bandwidth += rhs.bandwidth;
        utilization += rhs.utilization;
        latency += rhs.latency;
        return *this;
    }

    MonitoredSample& operator-=(const MonitoredSample& rhs){
        watts -= rhs.watts;
        bandwidth -= rhs.bandwidth;
        utilization -= rhs.utilization;
        latency -= rhs.latency;
        return *this;
    }

    MonitoredSample& operator*=(const MonitoredSample& rhs){
        watts *= rhs.watts;
        bandwidth *= rhs.bandwidth;
        utilization *= rhs.utilization;
        latency *= rhs.latency;
        return *this;
    }

    MonitoredSample& operator/=(const MonitoredSample& rhs){
        watts /= rhs.watts;
        bandwidth /= rhs.bandwidth;
        utilization /= rhs.utilization;
        latency /= rhs.latency;
        return *this;
    }

    MonitoredSample operator/=(double x){
        watts /= x;
        bandwidth /= x;
        utilization /= x;
        latency /= x;
        return *this;
    }

    MonitoredSample operator*=(double x){
        watts *= x;
        bandwidth *= x;
        utilization *= x;
        latency *= x;
        return *this;
    }
}MonitoredSample;

inline MonitoredSample operator+(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r += rhs;
    return r;
}

inline MonitoredSample operator-(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r -= rhs;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r *= rhs;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r /= rhs;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r /= x;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r *= x;
    return r;
}

inline ostream& operator<<(ostream& os, const MonitoredSample& obj){
    os << "[";
    os << "Watts: " << obj.watts << " ";
    os << "Bandwidth: " << obj.bandwidth << " ";
    os << "Utilization: " << obj.utilization << " ";
    os << "Latency: " << obj.latency << " ";
    os << "]";
    return os;
}

inline ofstream& operator<<(ofstream& os, const MonitoredSample& obj){
    os << obj.watts.cores << "\t";
    os << obj.bandwidth << "\t";
    os << obj.utilization << "\t";
    os << obj.latency << "\t";
    return os;
}

inline JoulesCpu squareRoot(const JoulesCpu& x){
    JoulesCpu r;
    r.cores = x.cores?sqrt(x.cores):0;
    r.cpu = x.cpu?sqrt(x.cpu):0;
    r.graphic = x.graphic?sqrt(x.graphic):0;
    r.dram = x.dram?sqrt(x.dram):0;
    return r;
}

inline MonitoredSample squareRoot(const MonitoredSample& x){
    MonitoredSample r;
    r.watts = squareRoot(x.watts);
    r.bandwidth = x.bandwidth?sqrt(x.bandwidth):0;
    r.utilization = x.utilization?sqrt(x.utilization):0;
    r.latency = x.latency?sqrt(x.latency):0;
    return r;
}


inline void regularize(MonitoredSample& x){
    if(x.watts.cpu < 0){x.watts.cpu = 0;}
    if(x.watts.cores < 0){x.watts.cores = 0;}
    if(x.watts.graphic < 0){x.watts.graphic = 0;}
    if(x.watts.dram < 0){x.watts.dram = 0;}
    if(x.bandwidth < 0){x.bandwidth = 0;}
    if(x.utilization < 0){x.utilization = 0;}
    if(x.latency < 0){x.latency = 0;}
}

}

#endif /* ADAPTIVE_FASTFLOW_FARM_HPP_ */
