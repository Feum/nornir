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

#include "knob.hpp"
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

//TODO REMOVE USING
using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

struct MonitoredSample;

typedef enum{
    KNOB_TYPE_WORKERS = 0,
    KNOB_TYPE_MAPPING,
    KNOB_TYPE_FREQUENCY,
    KNOB_TYPE_NUM // <---- This must always be the last value
}KnobType;

typedef enum{
    KNOB_VALUE_UNDEF = 0,
    KNOB_VALUE_REAL,
    KNOB_VALUE_RELATIVE
}KnobValueType;

class KnobsValues{
private:
    KnobValueType _type;
    double _values[KNOB_TYPE_NUM];
public:
    KnobsValues(KnobValueType type = KNOB_VALUE_UNDEF):_type(type){;}

    inline bool areRelative() const{return _type == KNOB_VALUE_RELATIVE;}

    inline bool areReal() const{return _type == KNOB_VALUE_REAL;}

    inline double& operator[](KnobType idx){
        return _values[idx];
    }

    inline double operator[](KnobType idx) const{
        return _values[idx];
    }
};

class FarmConfiguration: public mammut::utils::NonCopyable {
private:
    Knob* _knobs[KNOB_TYPE_NUM];
    const Parameters& _p;
    std::vector<KnobsValues> _combinations;
    void combinations(vector<vector<double> > array, size_t i,
                      vector<double> accum);
    void setRelativeValues(const KnobsValues& values);
    void setRealValues(const KnobsValues& values);
public:
    FarmConfiguration(const Parameters& p, ff::ff_farm<>& farm);

    ~FarmConfiguration();

    /**
     * Gets all the possible combinations of knobs values.
     * @return A vector containing all the possible combinations
     *         of knobs values.
     */
    const std::vector<KnobsValues>& getAllRealCombinations();

    /**
     * Sets the highest frequency to reduce the reconfiguration time.
     */
    void setFastReconfiguration();

    /**
     * Returns a specified knob.
     * @param t The type of the knob to return.
     * @return The specified knob.
     */
    const Knob* getKnob(KnobType t) const;

    /**
     * Sets all the knobs to their maximum.
     */
    void maxAllKnobs();

    /**
     * Returns the real value of a specific knob.
     * @param t The type of the knob.
     * @return The real value of the specified knob.
     */
    double getRealValue(KnobType t) const;

    /**
     * Returns the real values for all the knobs.
     * @return The real values for all the knobs.
     */
    KnobsValues getRealValues() const;

    /**
     * Sets values for the knobs (may be relative or real).
     * @param values The values of the knobs.
     */
    void setValues(const KnobsValues& values);

};


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

    // The vector of active workers.
    std::vector<AdaptiveNode*> _activeWorkers;

    // The current configuration of the farm.
    FarmConfiguration _configuration;

    // The voltage table.
    VoltageTable _voltageTable;

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

#ifdef DEBUG_MANAGER
    ofstream samplesFile;
#endif

    /**
     * Set a specified domain to the highest frequency.
     * @param domain The domain.
     */
    void setDomainToHighestFrequency(const Domain* domain);


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
     * Returns the voltage at a specific combinations of knobs values.
     * @param values The combination.
     * @return The voltage at that combination.
     */
    double getVoltage(const KnobsValues& values) const;

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
     * Computes the best relative knobs values for the farm.
     * @return The best relative knobs values.
     */
    KnobsValues getBestKnobsValues();

    /**
     * Checks if the application terminated.
     */
    bool terminated();

    /**
     * Changes the knobs.
     * @param values The new knobs values.
     */
    void changeKnobs(KnobsValues values);

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

    /**
     * Operations that need to take place before running the nodes.
     */
    void initNodesPreRun();

    /**
     * Operations that need to take place after running the nodes.
     */
    void initNodesPostRun();

    /**
     * Cleans the adaptive nodes.
     */
    void cleanNodes();
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

inline ostream& operator<<(ostream& os, const KnobsValues& obj){
    os << "[";
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        os << obj[(KnobType) i] << ", ";
    }
    os << "]";
    return os;
}

inline bool operator==(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(lhs[(KnobType) i] !=
           rhs[(KnobType) i]){
            return false;
        }
    }
    return true;
}

inline bool operator!=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator==(lhs,rhs);
}

inline bool operator<(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(lhs[(KnobType) i] <
           rhs[(KnobType) i]){
            return true;
        }else if(lhs[(KnobType) i] >
                 rhs[(KnobType) i]){
            return false;
        }
    }
    return false;
}

inline bool operator>(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    return operator< (rhs,lhs);
}

inline bool operator<=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator> (lhs,rhs);
}

inline bool operator>=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
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
