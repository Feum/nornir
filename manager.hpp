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

/*! \file farm.hpp
 * \brief Implementation of an adaptive fastflow farm.
 *
 * To let an existing fastflow farm-based adaptive, follow these steps:
 *  1. Emitter, Workers and Collector of the farm must extend
 *     adpff::adpff_node instead of ff::ff_node
 *  2. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyWorkersChange virtual
 *     method.
 */

#ifndef ADAPTIVE_FASTFLOW_FARM_HPP_
#define ADAPTIVE_FASTFLOW_FARM_HPP_

#define TRACE_FASTFLOW
#define FF_TASK_CALLBACK

#include "parameters.hpp"
#include "predictors.hpp"
#include "adpnode.hpp"
#include "utils.hpp"

#include <ff/farm.hpp>
#include <mammut/module.hpp>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

#include <cmath>
#include <iostream>
#include <limits>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_FARM
#define DEBUG(x) do { cerr << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace adpff{

class AdaptivityParameters;
class ManagerFarm;

using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

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
                                     uint durationMs){
        return ((double)cs.duration /
                (double)durationMs) * 100.0;
    }
public:
    Observer(string statsFile = "stats.csv",
             string calibrationFile = "calibration.csv",
             string summaryFile = "summary.csv"):
            _startMonitoringMs(0),
            _totalWatts(0),
            _totalBw(0),
            _numSamples(0){
        _statsFile.open(statsFile.c_str());
        _calibrationFile.open(calibrationFile.c_str());
        _summaryFile.open(summaryFile.c_str());
        if(!_statsFile.is_open() ||
           !_calibrationFile.is_open() ||
           !_summaryFile.is_open()){
            throw runtime_error("Observer: Impossible to open file.");
        }
        _statsFile << "TimestampMillisecs" << "\t";
        _statsFile << "[[EmitterVc][WorkersVc][CollectorVc]]" << "\t";
        _statsFile << "Workers" << "\t";
        _statsFile << "Frequency" << "\t";
        _statsFile << "CurrentBandwidth" << "\t";
        _statsFile << "SmoothedBandwidth" << "\t";
        _statsFile << "CoeffVarBandwidth" << "\t";
        _statsFile << "SmoothedUtilization" << "\t";
        _statsFile << "SmoothedWattsCpu" << "\t";
        _statsFile << "SmoothedWattsCores" << "\t";
        _statsFile << "SmoothedWattsGraphic" << "\t";
        _statsFile << "SmoothedWattsDram" << "\t";
        _statsFile << endl;

        _calibrationFile << "NumSteps" << "\t";
        _calibrationFile << "Duration" << "\t";
        _calibrationFile << "Time%" << "\t";
        _calibrationFile << endl;

        _summaryFile << "Watts" << "\t";
        _summaryFile << "Bandwidth" << "\t";
        _summaryFile << "CompletionTime" << "\t";
        _summaryFile << "Calibration%" << "\t";
        _summaryFile << endl;
    }

    virtual ~Observer(){
        _statsFile.close();
        _calibrationFile.close();
        _summaryFile.close();
    }

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
                         JoulesCpu smoothedWatts){
        _statsFile << timeStamp - _startMonitoringMs << "\t";
        _statsFile << "[";
        if(emitterVirtualCore){
            _statsFile << "[" << emitterVirtualCore->getVirtualCoreId() << "]";
        }

        _statsFile << "[";
        for(size_t i = 0; i < workersVirtualCore.size(); i++){
            _statsFile << workersVirtualCore.at(i)->getVirtualCoreId() << ",";
        }
        _statsFile << "]";

        if(collectorVirtualCore){
            _statsFile << "[" << collectorVirtualCore->getVirtualCoreId() << "]";
        }
        _statsFile << "]" << "\t";

        _statsFile << workers << "\t";
        _statsFile << frequency << "\t";
        _statsFile << currentBandwidth << "\t";
        _statsFile << smoothedBandwidth << "\t";
        _statsFile << coeffVarBandwidth << "\t";
        _statsFile << smoothedUtilization << "\t";

        _statsFile << smoothedWatts.cpu << "\t";
        _statsFile << smoothedWatts.cores << "\t";
        _statsFile << smoothedWatts.graphic << "\t";
        _statsFile << smoothedWatts.dram << "\t";

        _statsFile << endl;

        _totalWatts += smoothedWatts.cores;
        _totalBw += currentBandwidth;
        _numSamples++;
    }

    virtual void calibrationStats(const vector<CalibrationStats>&
                                  calibrationStats,
                                  uint durationMs){

        for(size_t i = 0; i < calibrationStats.size(); i++){
            const CalibrationStats& cs = calibrationStats.at(i);
            _calibrationFile << cs.numSteps << "\t";
            _calibrationFile << cs.duration << "\t";
            _calibrationFile << calibrationDurationToPerc(cs, durationMs) << "\t";
            _calibrationFile << endl;
        }
    }

    virtual void summaryStats(const vector<CalibrationStats>&
                              calibrationStats,
                              uint durationMs){
        double totalCalibrationPerc = 0.0;
        for(size_t i = 0; i < calibrationStats.size(); i++){
            const CalibrationStats& cs = calibrationStats.at(i);
            totalCalibrationPerc += calibrationDurationToPerc(cs, durationMs);
        }

        _summaryFile << _totalWatts / (double) _numSamples << "\t";
        _summaryFile << _totalBw / (double) _numSamples << "\t";
        _summaryFile << (double) durationMs / 1000.0 << "\t";
        _summaryFile << totalCalibrationPerc << "\t";
        _summaryFile << endl;
    }
};

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

    FarmConfiguration():numWorkers(0), frequency(0){;}
    FarmConfiguration(uint numWorkers, Frequency frequency = 0):
                     numWorkers(numWorkers), frequency(frequency){;}
}FarmConfiguration;

ostream& operator<<(ostream& os, const FarmConfiguration& obj){
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

ostream& operator<<(ostream& os, const MonitoredSample& obj){
    os << "[";
    os << "Watts: " << obj.watts << " ";
    os << "Bandwidth: " << obj.bandwidth << " ";
    os << "Utilization: " << obj.utilization << " ";
    os << "Latency: " << obj.latency << " ";
    os << "]";
    return os;
}

ofstream& operator<<(ofstream& os, const MonitoredSample& obj){
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

/*!
 * \internal
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
private:
    // The managed farm.
    ff_farm<>* _farm;

    // The parameters used to take management decisions.
    AdaptivityParameters _p;

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
    adpff_node* _emitter;

    // The collector (if present).
    adpff_node* _collector;

    // The currently running workers.
    vector<adpff_node*> _activeWorkers;

    // Workers that can run but are not currently running.
    vector<adpff_node*> _inactiveWorkers;

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
    bool reconfigureWorkers() const{
        return true;
    }

    /**
     * Returns true if the frequencies must be reconfigured.
     * @return true if the frequencies must be reconfigured,
     *         false otherwise
     */
    bool reconfigureFrequency() const{
        return _p.strategyFrequencies == STRATEGY_FREQUENCY_YES ||
               _p.strategyFrequencies == STRATEGY_FREQUENCY_MIN_CORES;
    }

    /**
     * Returns the number of dimensions of a configuration,
     * i.e. the number of different decisions that can
     * be taken at runtime.
     */
    uint getConfigurationDimension() const{
        uint numDimensions = 0;
        if(reconfigureWorkers()){++numDimensions;}
        if(reconfigureFrequency()){++numDimensions;}
        return numDimensions;
    }

    /**
     * If possible, finds a set of physical cores belonging to domains
     * different from those of virtual cores in 'virtualCores' vector.
     * @param virtualCores A vector of virtual cores.
     * @return A set of physical cores that can always run at the highest
     *         frequency.
     */
    vector<PhysicalCore*> getSepDomainsPhyCores(
            const vector<VirtualCore*>& virtualCores) const{
        vector<Domain*> allDomains = _cpufreq->getDomains();
        vector<Domain*> hypotheticWDomains = _cpufreq->getDomains(virtualCores);
        vector<PhysicalCore*> physicalCoresInUnusedDomains;
        if(allDomains.size() > hypotheticWDomains.size()){
           for(size_t i = 0; i < allDomains.size(); i++){
               Domain* currentDomain = allDomains.at(i);
               if(!contains(hypotheticWDomains, currentDomain)){
                   insertToEnd(_topology->virtualToPhysical(
                                      currentDomain->getVirtualCores()),
                                      physicalCoresInUnusedDomains);
               }
           }
        }
        return physicalCoresInUnusedDomains;
    }

    /**
     * Set a specified domain to the highest frequency.
     * @param domain The domain.
     */
    void setDomainToHighestFrequency(const Domain* domain){
        if(!domain->setGovernor(GOVERNOR_PERFORMANCE)){
            if(!domain->setGovernor(GOVERNOR_USERSPACE) ||
               !domain->setHighestFrequencyUserspace()){
                throw runtime_error("AdaptivityManagerFarm: Fatal error while "
                                    "setting highest frequency for sensitive "
                                    "emitter/collector. Try to run it without "
                                    "sensitivity parameters.");
            }
        }
    }

    /**
     * Computes the available virtual cores, sorting them according to
     * the specified mapping strategy.
     * @return The available virtual cores, sorted according to the
     *         specified mapping strategy.
     */
    vector<VirtualCore*> getAvailableVirtualCores(){
        vector<VirtualCore*> r;
        if(_p.strategyMapping == STRATEGY_MAPPING_AUTO){
            _p.strategyMapping = STRATEGY_MAPPING_LINEAR;
        }

        switch(_p.strategyMapping){
            case STRATEGY_MAPPING_LINEAR:{
               /*
                * Generates a vector of virtual cores to be used for linear
                * mapping.node. It contains first one virtual core per physical
                * core (virtual cores on the same CPU are consecutive).
                * Then, the other groups of virtual cores follow.
                */
                vector<Cpu*> cpus = _topology->getCpus();

                size_t virtualUsed = 0;
                size_t virtualPerPhysical;
                if(_p.strategyHyperthreading != STRATEGY_HT_NO){
                    virtualPerPhysical = _topology->getVirtualCores().size() /
                                         _topology->getPhysicalCores().size();
                }else{
                    virtualPerPhysical = 1;
                }
                while(virtualUsed < virtualPerPhysical){
                    for(size_t i = 0; i < cpus.size(); i++){
                        vector<PhysicalCore*> phyCores = cpus.at(i)->
                                                         getPhysicalCores();
                        for(size_t j = 0; j < phyCores.size(); j++){
                            r.push_back(phyCores.at(j)->getVirtualCores().
                                                        at(virtualUsed));
                        }
                    }
                    ++virtualUsed;
                }
            }break;
            case STRATEGY_MAPPING_CACHE_EFFICIENT:{
                throw runtime_error("Not yet supported.");
            }
            default:
                break;
        }
        return r;
    }

    /**
     * Returns the number of scalable service nodes.
     * @return The number of scalable service nodes.
     */
    uint numScalableServiceNodes(){
        return (_emitter &&
                _p.mappingEmitter != SERVICE_NODE_MAPPING_PERFORMANCE) +
               (_collector &&
                _p.mappingEmitter != SERVICE_NODE_MAPPING_PERFORMANCE);
    }

    /**
     * Manages mapping of emitter and collector.
     */
    void manageServiceNodesPerformance(){
        if((_p.mappingEmitter != SERVICE_NODE_MAPPING_ALONE) && !_emitter){
            _p.mappingEmitter = SERVICE_NODE_MAPPING_ALONE;
        }

        if((_p.mappingCollector != SERVICE_NODE_MAPPING_ALONE) && !_collector){
            _p.mappingCollector = SERVICE_NODE_MAPPING_ALONE;
        }

        if((_p.mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE &&
                !_emitterSensitivitySatisfied) ||
           (_p.mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE &&
                   !_collectorSensitivitySatisfied)){
            size_t scalableVCNum = _activeWorkers.size() +
                                   numScalableServiceNodes();

            // When sensitive is specified, we always choose the WEC mapping.
            vector<VirtualCore*>::const_iterator scalEnd;
            if(scalableVCNum < _availableVirtualCores.size()){
                scalEnd = _availableVirtualCores.begin() + scalableVCNum;
            }else{
                scalEnd = _availableVirtualCores.end();
            }

            vector<VirtualCore*> scalableVC(_availableVirtualCores.begin(),
                                            scalEnd);
            vector<PhysicalCore*> perfPhyCores;
            perfPhyCores = getSepDomainsPhyCores(scalableVC);
            if(perfPhyCores.size()){
                size_t index = 0;

                if(_p.mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE){
                    VirtualCore* vc = perfPhyCores.at(index)->getVirtualCore();
                    setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                    _emitterVirtualCore = vc;
                    _emitterSensitivitySatisfied = true;
                    index = (index + 1) % perfPhyCores.size();
                }

                if(_p.mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE){
                    VirtualCore* vc = perfPhyCores.at(index)->getVirtualCore();
                    setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                    _collectorVirtualCore = vc;
                    _collectorSensitivitySatisfied = true;
                    index = (index + 1) % perfPhyCores.size();
                }
            }
        }
    }

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
                           size_t& collectorIndex){
        size_t nextIndex = 0;
        if(_emitter && !_emitterVirtualCore){
            emitterIndex = nextIndex;
            if(_p.mappingEmitter != SERVICE_NODE_MAPPING_COLLAPSED){
                nextIndex = (nextIndex + 1) % _availableVirtualCores.size();
            }
        }

        firstWorkerIndex = nextIndex;
        nextIndex = (nextIndex + _activeWorkers.size()) %
                    _availableVirtualCores.size();

        if(_collector && !_collectorVirtualCore){
            if(_p.mappingCollector == SERVICE_NODE_MAPPING_COLLAPSED){
                nextIndex = (nextIndex - 1) % _availableVirtualCores.size();
            }
            collectorIndex = nextIndex;
        }
    }

    /**
     * Computes the virtual cores where the nodes must be mapped
     * and pins these nodes on the virtual cores.
     */
    void mapNodesToVirtualCores(){
        size_t emitterIndex, firstWorkerIndex, collectorIndex;
        getMappingIndexes(emitterIndex, firstWorkerIndex, collectorIndex);

        //TODO: Che succede se la farm ha l'emitter di default? (non estende
        //      adaptive node e quindi la getThreadHandler non c'e' (fallisce
        //      lo static_cast prima)):(
        if(_emitter){
            if(!_emitterVirtualCore){
                _emitterVirtualCore = _availableVirtualCores.at(emitterIndex);
            }
            if(!_emitterVirtualCore->isHotPlugged()){
                _emitterVirtualCore->hotPlug();
            }
            _emitter->move(_emitterVirtualCore);
        }

        for(size_t i = 0; i < _activeWorkers.size(); i++){
            VirtualCore* vc = _availableVirtualCores.at((firstWorkerIndex + i) %
                              _availableVirtualCores.size());
            _activeWorkersVirtualCores.push_back(vc);
            if(!vc->isHotPlugged()){
                vc->hotPlug();
            }
            _activeWorkers.at(i)->move(vc);
        }

        if(_collector){
            if(!_collectorVirtualCore){
                _collectorVirtualCore = _availableVirtualCores.
                                        at(collectorIndex);
            }
            if(!_collectorVirtualCore->isHotPlugged()){
                _collectorVirtualCore->hotPlug();
            }
            _collector->move(_collectorVirtualCore);
        }
    }

    /**
     * Applies the OFF strategy for a specified set of virtual cores.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategyOff(const vector<VirtualCore*>& unusedVc){
        for(size_t i = 0; i < unusedVc.size(); i++){
            VirtualCore* vc = unusedVc.at(i);
            if(vc->isHotPluggable() && vc->isHotPlugged()){
                vc->hotUnplug();
            }
        }
    }

    /**
     * Applies the LOWEST_FREQUENCY strategy for a specified set of
     * virtual cores.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategyLowestFreq(const vector<VirtualCore*>& vc){
        vector<Domain*> unusedDomains = _cpufreq->getDomainsComplete(vc);
        for(size_t i = 0; i < unusedDomains.size(); i++){
            Domain* domain = unusedDomains.at(i);
            if(!domain->setGovernor(GOVERNOR_POWERSAVE)){
                if(!domain->setGovernor(GOVERNOR_USERSPACE) ||
                   !domain->setLowestFrequencyUserspace()){
                    throw runtime_error("AdaptivityManagerFarm: Impossible to "
                                        "set lowest frequency for unused "
                                        "virtual cores.");
                }
            }
        }
    }

    /**
     * Applies a specific strategy for a specified set of virtual cores.
     * @param strategyUnused The unused strategy.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVCStrategy(StrategyUnusedVirtualCores strategyUnused,
                               const vector<VirtualCore*>& vc){
        switch(strategyUnused){
            case STRATEGY_UNUSED_VC_OFF:{
                applyUnusedVCStrategyOff(vc);
            }break;
            case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                applyUnusedVCStrategyLowestFreq(vc);
            }break;
            default:{
                return;
            }
        }
    }

    /**
     * Apply the strategies for inactive and unused virtual cores.
     */
    void applyUnusedVCStrategy(){
        /**
         * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the
         * virtual cores on a domain, we can also lower its frequency to
         * the minimum.
         */
        vector<VirtualCore*> virtualCores;
        if(_p.strategyInactiveVirtualCores != STRATEGY_UNUSED_VC_NONE){
            insertToEnd(_inactiveWorkersVirtualCores, virtualCores);
        }
        if(_p.strategyUnusedVirtualCores != STRATEGY_UNUSED_VC_NONE){
            insertToEnd(_unusedVirtualCores, virtualCores);
        }
        applyUnusedVCStrategy(STRATEGY_UNUSED_VC_LOWEST_FREQUENCY, virtualCores);


        virtualCores.clear();
        if(_p.strategyInactiveVirtualCores == STRATEGY_UNUSED_VC_OFF){
            insertToEnd(_inactiveWorkersVirtualCores, virtualCores);
        }
        if(_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_OFF){
            insertToEnd(_unusedVirtualCores, virtualCores);
        }
        applyUnusedVCStrategy(STRATEGY_UNUSED_VC_OFF, virtualCores);
    }

    /**
     * Updates the scalable domains vector.
     */
    void updateScalableDomains(){
        vector<VirtualCore*> scalableVirtualCores = _activeWorkersVirtualCores;
        /**
         * Node sensitivity may be not satisfied both because it was not
         * requested, or because it was requested but it was not possible
         * to satisfy it. In both cases, we need to scale the virtual core
         * of the node as we scale the others.
         */
        if(_emitter && !_emitterSensitivitySatisfied){
            scalableVirtualCores.push_back(_emitterVirtualCore);
        }
        if(_collector && !_collectorSensitivitySatisfied){
            scalableVirtualCores.push_back(_collectorVirtualCore);
        }

        _scalableDomains = _cpufreq->getDomains(scalableVirtualCores);
    }

    /**
     * Set a specific P-state for the virtual cores used by
     * the current active workers, emitter (if not sensitive) and
     * collector (if not sensitive).
     * @param frequency The frequency to be set.
     */
    void updatePstate(Frequency frequency){
        updateScalableDomains();
        Domain* currentDomain;
        for(size_t i = 0; i < _scalableDomains.size(); i++){
            currentDomain = _scalableDomains.at(i);
            if(!currentDomain->setGovernor(_p.frequencyGovernor)){
                throw runtime_error("AdaptivityManagerFarm: Impossible to "
                                         "set the specified governor.");
            }
            if(_p.frequencyGovernor != GOVERNOR_USERSPACE){
                if(!currentDomain->setGovernorBounds(_p.frequencyLowerBound,
                                                     _p.frequencyUpperBound)){
                    throw runtime_error("AdaptivityManagerFarm: Impossible "
                                             "to set the specified governor's "
                                             "bounds.");
                }
            }else if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                if(!currentDomain->setFrequencyUserspace(frequency)){
                    throw runtime_error("AdaptivityManagerFarm: Impossible "
                                             "to set the specified frequency.");
                }
            }
        }
    }

    /**
     * Map the nodes to virtual cores and
     * prepares frequencies and governors for running.
     */
    void mapAndSetFrequencies(){
        if(_p.strategyMapping == STRATEGY_MAPPING_NO){
            return;
        }

        manageServiceNodesPerformance();
        mapNodesToVirtualCores();
        for(size_t i = 0; i < _availableVirtualCores.size(); i++){
            VirtualCore* vc = _availableVirtualCores.at(i);
            if(vc != _emitterVirtualCore && vc != _collectorVirtualCore &&
               !contains(_activeWorkersVirtualCores, vc) &&
               !contains(_inactiveWorkersVirtualCores, vc)){
                _unusedVirtualCores.push_back(vc);
            }
        }
        updateUsedCpus();
        applyUnusedVCStrategy();

        // Insert dummy constant frequency
        _availableFrequencies.push_back(1.0);

        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO){
            if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                // We suppose that all the domains have the same
                // available frequencies.
                _availableFrequencies = _cpufreq->getDomains().at(0)->
                                        getAvailableFrequencies();

                // Remove turbo boost frequency.
                if(!_p.turboBoost){
                    if(intToString(_availableFrequencies.back()).at(3) == '1'){
                        _availableFrequencies.pop_back();
                    }
                }

                // Sets the current frequency to the highest possible.
                _currentConfiguration.frequency = _availableFrequencies.back();
            }
            updatePstate(_currentConfiguration.frequency);
        }
    }

    /**
     * Returns the maximum primary prediction error.
     * @return The maximum primary prediction error.
     */
    double getMaxPredictionErrorPrimary() const{
        double r = _p.maxPrimaryPredictionError;
        if(_p.strategyPredictionErrorPrimary ==
           STRATEGY_PREDICTION_ERROR_COEFFVAR){
            r = max(r, getPrimaryValue(_samples->coefficientVariation()));
        }
        return r;
    }

    /**
     * Returns the maximum secondary prediction error.
     * @return The maximum secondary prediction error.
     */
    double getMaxPredictionErrorSecondary() const{
        double r = _p.maxSecondaryPredictionError;
        if(_p.strategyPredictionErrorSecondary ==
           STRATEGY_PREDICTION_ERROR_COEFFVAR){
            r = max(r, getSecondaryValue(_samples->coefficientVariation()));
        }
        return r;
    }

    /**
     * Returns the primary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The primary value of a sample according to
     * the required contract.
     */
    double getPrimaryValue(const MonitoredSample& sample) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return sample.utilization;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return sample.bandwidth;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return sample.watts.cores;
            }break;
            default:{
                return 0;
            }break;
        }
    }

    /**
     * Returns the secondary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The secondary value of a sample according to
     * the required contract.
     */
    double getSecondaryValue(const MonitoredSample& sample) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return sample.watts.cores;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return sample.bandwidth;
            }break;
            default:{
                return 0;
            }break;
        }
    }

    /**
     * Returns the primary value according to the required contract.
     * @return The primary value according to the required contract.
     */
    double getPrimaryValue() const{
        return getPrimaryValue(_samples->average());
    }

    /**
     * Returns the secondary value according to the required contract.
     * @return The secondary value according to the required contract.
     */
    double getSecondaryValue() const{
        return getSecondaryValue(_samples->average());
    }

    /**
     * Checks if the contract is violated.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated() const{
        double tolerance = 0;
        double maxError = getMaxPredictionErrorPrimary();
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                tolerance = ((_p.overloadThresholdFarm -
                              _p.underloadThresholdFarm) * maxError) / 100.0;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                tolerance = (_p.requiredBandwidth * maxError) / 100.0;
            }break;
            case CONTRACT_POWER_BUDGET:{
                tolerance = (_p.powerBudget * maxError) / 100.0;
            }break;
            default:{
                return false;
            }break;
        }
        return !isFeasiblePrimaryValue(getPrimaryValue(), tolerance);
    }

    bool isFeasiblePrimaryValue(double value, double tolerance = 0) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return value > _p.underloadThresholdFarm - tolerance &&
                       value < _p.overloadThresholdFarm + tolerance;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return value > _p.requiredBandwidth - tolerance;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return value < _p.powerBudget + tolerance;
            }break;
            default:{
                return false;
            }break;
        }
        return false;
    }

    /**
     * Returns the voltage at a specific configuration.
     * @param configuration The configuration.
     * @return The voltage at that configuration.
     */
    double getVoltage(const FarmConfiguration& configuration) const{
        VoltageTableKey key(configuration.numWorkers,
                                     configuration.frequency);
        VoltageTableIterator it = _voltageTable.find(key);
        if(it != _voltageTable.end()){
            return it->second;
        }else{
            throw runtime_error("Frequency and/or number of virtual cores "
                                     "not found in voltage table.");
        }
    }

    /**
     * Checks if x is a best suboptimal monitored value than y.
     * @param x The first monitored value.
     * @param y The second monitored value.
     * @return True if x is a best suboptimal monitored value than y,
     *         false otherwise.
     */
    bool isBestSuboptimalValue(double x, double y) const{
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

    /**
     * Returns true if x is a best secondary value than y, false otherwise.
     */
    bool isBestSecondaryValue(double x, double y) const{
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

    /**
     * Computes the new configuration of the farm after a contract violation.
     * @return The new configuration.
     */
    FarmConfiguration getNewConfiguration(){
        FarmConfiguration bestConfiguration;
        FarmConfiguration bestSuboptimalConfiguration = _currentConfiguration;

        double primaryPrediction = 0;
        double secondaryPrediction = 0;

        double bestPrimaryPrediction = 0;
        double bestSecondaryPrediction = 0;
        double primaryPredictionSub = 0;
        double secondaryPredictionSub = 0;

        double bestSecondaryValue = 0;
        double bestSuboptimalValue = getPrimaryValue();

        bool feasibleSolutionFound = false;

        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                // We have to minimize the power/energy.
                bestSecondaryValue = numeric_limits<double>::max();
            }break;
            case CONTRACT_POWER_BUDGET:{
                // We have to maximize the bandwidth.
                bestSecondaryValue = numeric_limits<double>::min();
            }break;
            default:{
                ;
            }break;
        }

        _primaryPredictor->prepareForPredictions();
        _secondaryPredictor->prepareForPredictions();

        unsigned int remainingTime = 0;
        for(size_t i = 1; i <= _maxNumWorkers; i++){
            for(size_t j = 0; j < _availableFrequencies.size(); j++){
                FarmConfiguration currentConf(i, _availableFrequencies.at(j));
                primaryPrediction = _primaryPredictor->predict(currentConf);
                switch(_p.contractType){
                    case CONTRACT_PERF_COMPLETION_TIME:{
                        remainingTime = (double) _remainingTasks /
                                         primaryPrediction;
                    }break;
                    case CONTRACT_PERF_UTILIZATION:{
                        primaryPrediction = (_samples->average().bandwidth /
                                             primaryPrediction) *
                                             _samples->average().utilization;
                    }break;
                    default:{
                        ;
                    }
                }

                if(isFeasiblePrimaryValue(primaryPrediction)){
                    secondaryPrediction = _secondaryPredictor->
                                          predict(currentConf);
                    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                        secondaryPrediction *= remainingTime;
                    }
                    if(isBestSecondaryValue(secondaryPrediction,
                                            bestSecondaryValue)){
                        bestSecondaryValue = secondaryPrediction;
                        bestConfiguration = currentConf;
                        feasibleSolutionFound = true;
                        bestPrimaryPrediction = primaryPrediction;
                        bestSecondaryPrediction = secondaryPrediction;
                    }
                }else if(!feasibleSolutionFound &&
                         isBestSuboptimalValue(primaryPrediction,
                                               bestSuboptimalValue)){
                    bestSuboptimalValue = primaryPrediction;
                    bestSuboptimalConfiguration = currentConf;
                    primaryPredictionSub = primaryPrediction;
                    secondaryPredictionSub = secondaryPrediction;
                }
            }
        }

        if(feasibleSolutionFound){
            _primaryPrediction = bestPrimaryPrediction;
            _secondaryPrediction = bestSecondaryPrediction;
            return bestConfiguration;
        }else{
            _primaryPrediction = primaryPredictionSub;
            _secondaryPrediction = secondaryPredictionSub;
            return bestSuboptimalConfiguration;
        }
    }

    /**
     * Updates the currently used and unused CPUs.
     */
    void updateUsedCpus(){
        _usedCpus.clear();
        _unusedCpus.clear();
        for(size_t i = 0; i < _activeWorkersVirtualCores.size(); i++){
            CpuId cpuId = _activeWorkersVirtualCores.at(i)->getCpuId();
            if(!contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }
        if(_emitterVirtualCore){
            CpuId cpuId = _emitterVirtualCore->getCpuId();
            if(!contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }
        if(_collectorVirtualCore){
            CpuId cpuId = _collectorVirtualCore->getCpuId();
            if(!contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }

        vector<Cpu*> cpus = _topology->getCpus();
        for(size_t i = 0; i < cpus.size(); i++){
            if(!contains(_usedCpus, cpus.at(i)->getCpuId())){
                _unusedCpus.push_back(cpus.at(i)->getCpuId());
            }
        }
    }

    /**
     * Changes the active and inactive nodes according to the new configuration.
     * @param configuration The new configuration.
     */
    void changeActiveNodes(FarmConfiguration& configuration){
        uint workersNumDiff = abs(_currentConfiguration.numWorkers -
                                  configuration.numWorkers);
        if(_currentConfiguration.numWorkers > configuration.numWorkers){
            /** Move workers from active to inactive. **/
            moveEndToFront(_activeWorkers, _inactiveWorkers, workersNumDiff);
            moveEndToFront(_activeWorkersVirtualCores,
                           _inactiveWorkersVirtualCores, workersNumDiff);
        }else{
            /** Move workers from inactive to active. **/
            /**
             * We need to map them again because if virtual cores were
             * shutdown, the threads that were running on them have been moved
             * on different virtual cores. Accordingly, when started again they
             * would not be run on the correct virtual cores.
             */
            for(size_t i = 0; i < workersNumDiff; i++){
                VirtualCore* vc = _inactiveWorkersVirtualCores.at(i);
                if(!vc->isHotPlugged()){
                    vc->hotPlug();
                }
                _inactiveWorkers.at(i)->move(vc);
            }
            moveFrontToEnd(_inactiveWorkers, _activeWorkers, workersNumDiff);
            moveFrontToEnd(_inactiveWorkersVirtualCores,
                           _activeWorkersVirtualCores, workersNumDiff);
        }
    }

    /**
     * Changes the current farm configuration.
     * @param configuration The new configuration.
     * @param refine True if the data on the last configuration must be used
     *               to refine the model.
     */
    void changeConfiguration(FarmConfiguration configuration,
                             bool refine = true){
        if(configuration == _currentConfiguration){
            return;
        }

        if(!configuration.numWorkers){
            throw runtime_error("AdaptivityManagerFarm: fatal error, trying to "
                                "activate zero workers.");
        }else if(configuration.numWorkers > _maxNumWorkers){
            throw runtime_error("AdaptivityManagerFarm: fatal error, trying to "
                                "activate more workers than the maximum "
                                "allowed.");
        }

        /****************** Refine the model ******************/
        if(refine){
            _primaryPredictor->refine();
            _secondaryPredictor->refine();
        }

        /****************** Workers change started ******************/
        if(_currentConfiguration.numWorkers != configuration.numWorkers){
            vector<RollbackPoint> rollbackPoints;
            if(_p.fastReconfiguration){
                for(size_t i = 0; i < _scalableDomains.size(); i++){
                    Domain* d = _scalableDomains.at(i);
                    rollbackPoints.push_back(d->getRollbackPoint());
                    setDomainToHighestFrequency(d);
                }
            }

            /** Freezes farm. **/
            DEBUG("Asking the farm to freeze");
            _emitter->freezeAll((void*)(_collector?FF_EOSW:FF_GO_OUT));
            DEBUG("Wait for freezing");
            for(size_t i = 0; i < _activeWorkers.size(); i++){
                _activeWorkers.at(i)->wait_freezing();
            }
            if(_collector){
                _collector->wait_freezing();
            }
            DEBUG("Farm freezed");

            changeActiveNodes(configuration);
            updateUsedCpus();

            /** Notify the nodes that a reconfiguration is happening. **/
            if(_emitter){
                _emitter->notifyWorkersChange(_currentConfiguration.numWorkers,
                                              configuration.numWorkers);
            }
            for(size_t i = 0; i < configuration.numWorkers; i++){
                adpff_node* w = _activeWorkers.at(i);
                w->notifyWorkersChange(_currentConfiguration.numWorkers,
                                       configuration.numWorkers);
            }
            if(_collector){
                _collector->notifyWorkersChange(_currentConfiguration.numWorkers,
                                                configuration.numWorkers);
            }

            /** Prepare all the nodes for the new run. **/
            if(_emitter){
                _emitter->prepareToRun();
            }
            for(size_t i = 0; i < configuration.numWorkers; i++){
                _activeWorkers.at(i)->prepareToRun();
            }
            if(_collector){
                _collector->prepareToRun();
            }


            /** Start the farm again. **/
            DEBUG("Asking the farm to start again");
            _emitter->thawAll(configuration.numWorkers);
            if(_collector){
                _collector->thaw(true, configuration.numWorkers);
            }
            DEBUG("Farm started");
            if(_p.fastReconfiguration){
                _cpufreq->rollback(rollbackPoints);
            }
        }
        /****************** Workers change terminated ******************/
        applyUnusedVCStrategy();
        /****************** P-state change started ******************/
        //TODO: Maybe sensitivity could not be satisfied with the maximum
        // number of workers but could be satisfied now.
        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO){
            updatePstate(configuration.frequency);
        }
        /****************** P-state change terminated ******************/
        _currentConfiguration = configuration;

        /****************** Clean state ******************/
        _lastStoredSampleMs = getMillisecondsTime();
        _samples->reset();
        _energy->resetCountersCpu();
        _totalTasks = 0;
    }

    /**
     * Send data to observer.
     **/
    void observe(){
        if(_p.observer){
            _p.observer->observe(_lastStoredSampleMs,
                                 _currentConfiguration.numWorkers,
                                 _currentConfiguration.frequency,
                                 _emitterVirtualCore,
                                 _activeWorkersVirtualCores,
                                 _collectorVirtualCore,
                                 _samples->getLastSample().bandwidth,
                                 _samples->average().bandwidth,
                                 _samples->coefficientVariation().bandwidth,
                                 _samples->average().utilization,
                                 _samples->average().watts);
        }
    }

    /**
     * Asks the workers for their samples.
     */
    void askForWorkersSamples(){
        for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
            _activeWorkers.at(i)->askForSample();
        }
    }

    /**
     * Obtain workers samples.
     * @param sample A worker sample. It will be filled by this call with the
     *               global data of the farm.
     * @return True if all the workers are still running, false otherwise.
     */
    bool getWorkersSamples(WorkerSample& sample){
        adpff_node* w;
        sample = WorkerSample();
        for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
            WorkerSample tmp;
            w = _activeWorkers.at(i);
            if(!w->getSampleResponse(tmp,
                                     _p.strategyPolling,
                                     _samples->average().latency)){
                return false;
            }
            sample += tmp;
        }
        sample.loadPercentage /= _currentConfiguration.numWorkers;
        sample.latency /= _currentConfiguration.numWorkers;
        return true;
    }

    /**
     * Store a new sample.
     * @param askWorkers If false, the request to the worker is not sent.
     * @return false if the farm is not running anymore, true otherwise.
     **/
    bool storeNewSample(){
        askForWorkersSamples();

        MonitoredSample sample;
        WorkerSample ws;
        JoulesCpu joules;
        if(!getWorkersSamples(ws)){
            return false;
        }

        _totalTasks += ws.tasksCount;
        if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
            if(_remainingTasks > ws.tasksCount){
                _remainingTasks -= ws.tasksCount;
            }else{
                _remainingTasks = 0;
            }
        }

        for(size_t i = 0; i < _usedCpus.size(); i++){
            CounterCpu* cc = _energy->getCounterCpu(_usedCpus.at(i));
            joules += cc->getJoules();
        }
        for(size_t i = 0; i < _unusedCpus.size(); i++){
            CounterCpu* cc = _energy->getCounterCpu(_unusedCpus.at(i));
            joules += cc->getJoules();
        }

        double now = getMillisecondsTime();
        double durationSecs = (now - _lastStoredSampleMs) / 1000.0;
        _lastStoredSampleMs = now;

        sample.watts = joules / durationSecs;
        sample.utilization = ws.loadPercentage;
        // ATTENTION: Bandwidth is not the number of task since the
        //            last observation but the number of expected
        //            tasks that will be processed in 1 second.
        //            For this reason, if we sum all the bandwidths in
        //            the result observation file, we may have an higher
        //            number than the number of tasks.
        sample.bandwidth = ws.bandwidthTotal;
        sample.latency = ws.latency;

        _energy->resetCountersCpu();
        _samples->add(sample);

        DEBUGB(samplesFile << *_samples << "\n");
        return true;
    }

    /**
     * Returns true if the manager doesn't have still to check for a new
     * configuration.
     * @return True if the manager doesn't have still to check for a new
     * configuration.
     */
    bool persist() const{
        bool r = false;
        switch(_p.strategyPersistence){
            case STRATEGY_PERSISTENCE_SAMPLES:{
                r = _samples->size() < _p.persistenceValue;
            }break;
            case STRATEGY_PERSISTENCE_TASKS:{
                r = _totalTasks < _p.persistenceValue;
            }break;
            case STRATEGY_PERSISTENCE_VARIATION:{
                const MonitoredSample& variation =
                        _samples->coefficientVariation();
                r = getPrimaryValue(variation) < _p.persistenceValue &&
                    getSecondaryValue(variation) < _p.persistenceValue;
            }break;
        }
        return r;
    }


    /**
     * Initializes the calibrator.
     */
    void initCalibrator(){
        if(_p.strategyCalibration == STRATEGY_CALIBRATION_RANDOM){
            ;
        }else{
            _calibrator = new CalibratorLowDiscrepancy(*this);
        }
    }

    /**
     * Initializes the predictors
     */
    void initPredictors(){
        PredictorType primary, secondary;
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                primary = PREDICTION_BANDWIDTH;
                secondary = PREDICTION_POWER;
            }break;
            case CONTRACT_POWER_BUDGET:{
                primary = PREDICTION_POWER;
                secondary = PREDICTION_BANDWIDTH;
            }break;
            default:{
                return;
            }break;
        }

        switch(_p.strategyPrediction){
            case STRATEGY_PREDICTION_SIMPLE:{
                _primaryPredictor = new PredictorSimple(primary, *this);
                _secondaryPredictor = new PredictorSimple(secondary, *this);
                _calibrator = NULL;
            }break;
            case STRATEGY_PREDICTION_REGRESSION_LINEAR:{
                _primaryPredictor = new PredictorLinearRegression(primary,
                                                                  *this);
                _secondaryPredictor = new PredictorLinearRegression(secondary,
                                                                    *this);
                initCalibrator();
            }break;
            default:{
                ;
            }break;
        }
    }

public:
    /**
     * Creates a farm adaptivity manager.
     * @param farm The farm to be managed.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFarm(ff_farm<>* farm, AdaptivityParameters adaptivityParameters):
            _farm(farm),
            _p(adaptivityParameters),
            _startTimeMs(0),
            _cpufreq(_p.mammut.getInstanceCpuFreq()),
            _energy(_p.mammut.getInstanceEnergy()),
            _task(_p.mammut.getInstanceTask()),
            _topology(_p.mammut.getInstanceTopology()),
            _numCpus(_topology->getCpus().size()),
            _numPhysicalCores(_topology->getPhysicalCores().size()),
            _numPhysicalCoresPerCpu(_topology->getCpu(0)->getPhysicalCores().
                                    size()),
            _numVirtualCoresPerPhysicalCore(_topology->getPhysicalCore(0)->
                                            getVirtualCores().size()),
            _emitterSensitivitySatisfied(false),
            _collectorSensitivitySatisfied(false),
            _availableVirtualCores(getAvailableVirtualCores()),
            _emitterVirtualCore(NULL),
            _collectorVirtualCore(NULL){
        /** If voltage table file is specified, then load the table. **/
        if(_p.archData.voltageTableFile.compare("")){
            loadVoltageTable(_voltageTable,
                                      _p.archData.voltageTableFile);
        }

        _samples = NULL;
        switch(_p.strategySmoothing){
        case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
            _samples = new MovingAverageSimple<MonitoredSample>
                                    (_p.smoothingFactor);
        }break;
        case STRATEGY_SMOOTHING_EXPONENTIAL:{
            _samples = new MovingAverageExponential<MonitoredSample>
                                    (_p.smoothingFactor);
        }break;
        }

        DEBUGB(samplesFile.open("samples.csv"));
    }

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFarm(){
        delete _samples;
        if(_primaryPredictor){
            delete _primaryPredictor;
        }
        if(_secondaryPredictor){
            delete _secondaryPredictor;
        }
        if(_calibrator){
            delete _calibrator;
        }
        DEBUGB(samplesFile.close());
    }

    /**
     * Function executed by this thread.
     */
    void run(){
        _emitter = static_cast<adpff_node*>(_farm->getEmitter());
        _collector = static_cast<adpff_node*>(_farm->getCollector());
        svector<ff_node*> w = _farm->getWorkers();
        for(size_t i = 0; i < w.size(); i++){
            _activeWorkers.push_back(static_cast<adpff_node*>(w[i]));
        }
        _maxNumWorkers = _activeWorkers.size();
        _currentConfiguration = FarmConfiguration(_maxNumWorkers, 0);

        _farm->run_then_freeze(_maxNumWorkers);

        //TODO: Is already running when exits from that call?

        for(size_t i = 0; i < _activeWorkers.size(); i++){
            _activeWorkers.at(i)->init(_p.mammut, _p.archData.ticksPerNs);
        }
        if(_emitter){
            _emitter->init(_p.mammut, _p.archData.ticksPerNs);
        }
        if(_collector){
            _collector->init(_p.mammut, _p.archData.ticksPerNs);
        }

        _startTimeMs = getMillisecondsTime();

        if(_cpufreq->isBoostingSupported()){
            if(_p.turboBoost){
                _cpufreq->enableBoosting();
            }else{
                _cpufreq->disableBoosting();
            }
        }

        mapAndSetFrequencies();
        _energy->resetCountersCpu();

        _lastStoredSampleMs = _startTimeMs;
        if(_p.observer){
            _p.observer->_startMonitoringMs = _lastStoredSampleMs;
        }

        if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
            _remainingTasks = _p.expectedTasksNumber;
            _deadline = getMillisecondsTime()/1000.0 +
                        _p.requiredCompletionTime;
        }

        initPredictors();

        double microsecsSleep = 0;
        if(_p.contractType == CONTRACT_NONE){
            _farm->wait();
            storeNewSample();
            observe();
        }else{
            /* Force the first calibration point. **/
            if(_calibrator){
                changeConfiguration(_calibrator->getNextConfiguration(), false);
            }

            double startSample = getMillisecondsTime();
            while(true){
                double overheadMs = getMillisecondsTime() - startSample;
                microsecsSleep = ((double)_p.samplingInterval - overheadMs)*
                                  (double)MAMMUT_MICROSECS_IN_MILLISEC;
                if(microsecsSleep < 0){
                    microsecsSleep = 0;
                }
                usleep(microsecsSleep);
                startSample = getMillisecondsTime();

                if(!storeNewSample()){
                    goto controlLoopEnd;
                }

                if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                    uint now = getMillisecondsTime()/1000.0;
                    if(now >= _deadline){
                        _p.requiredBandwidth = numeric_limits<double>::max();
                    }else{
                        _p.requiredBandwidth = _remainingTasks /
                                               (_deadline - now);
                    }
                }

                observe();

                if(!persist()){
                    bool reconfigurationRequired = false;
                    FarmConfiguration nextConfiguration;

                    if(_calibrator){
                        nextConfiguration = _calibrator->getNextConfiguration();
                        if(nextConfiguration != _currentConfiguration){
                            reconfigurationRequired = true;
                        }
                    }else if(isContractViolated()){
                        nextConfiguration = getNewConfiguration();
                        reconfigurationRequired = true;
                    }

                    if(reconfigurationRequired){
                        changeConfiguration(nextConfiguration);
                        startSample = getMillisecondsTime();
                    }
                }
            }
        }
    controlLoopEnd:

        uint duration = getMillisecondsTime() - _startTimeMs;
        if(_p.observer){
            vector<CalibrationStats> cs;
            if(_calibrator){
                cs = _calibrator->getCalibrationsStats();
                _p.observer->calibrationStats(cs, duration);
            }
            _p.observer->summaryStats(cs, duration);
        }
    }
};

}

#endif /* ADAPTIVE_FASTFLOW_FARM_HPP_ */
