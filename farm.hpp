/*
 * farm.hpp
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
 *     adpff::adp_ff_node instead of ff::ff_node
 *  2. Replace the following calls (if present) in the farm nodes:
 *          svc           -> adp_svc
 *          svc_init      -> adp_svc_init
 *          svc_end       -> adp_svc_end
 *  3. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyWorkersChange virtual
 *     method.
 *  4. Substitute ff::ff_farm with adpff::adp_ff_farm. The maximum number of
 *     workers that can be activated correspond to the number of workers
 *     specified during farm creation.
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

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_FARM
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace adpff{

class AdaptivityParameters;
class AdaptivityManagerFarm;

using namespace ff;
using namespace mammut;

/*!
 * This class can be used to obtain statistics about reconfigurations
 * performed by the manager.
 * It can be extended by a user defined class to customize action to take
 * for each observed statistic.
 */
class Observer{
    friend class AdaptivityManagerFarm;
private:
    std::ofstream _statsFile;
    std::ofstream _calibrationFile;
    unsigned int _startMonitoringMs;
public:
    Observer(std::string statsFile = "stats.csv",
             std::string calibrationFile = "calibration.csv"):
            _startMonitoringMs(0){
        _statsFile.open(statsFile.c_str());
        _calibrationFile.open(calibrationFile.c_str());
        if(!_statsFile.is_open()){
            throw std::runtime_error("Observer: Impossible to open file.");
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
        _statsFile << std::endl;

        _calibrationFile << "NumSteps" << "\t";
        _calibrationFile << "Duration" << "\t";
        _calibrationFile << "Time%" << "\t";
        _calibrationFile << std::endl;
    }

    virtual ~Observer(){
        _statsFile.close();
        _calibrationFile.close();
    }

    virtual void observe(unsigned int timeStamp,
                         size_t workers,
                         cpufreq::Frequency frequency,
                         const topology::VirtualCore* emitterVirtualCore,
                         const std::vector<topology::VirtualCore*>& workersVirtualCore,
                         const topology::VirtualCore* collectorVirtualCore,
                         double currentBandwidth,
                         double smoothedBandwidth,
                         double coeffVarBandwidth,
                         double smoothedUtilization,
                         energy::JoulesCpu smoothedWatts){
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

        _statsFile << std::endl;
    }

    virtual void calibrationStats(const std::vector<CalibrationStats>&
                                  calibrationStats,
                                  uint totalDurationMs){
        double timePerc = 0.0;
        for(size_t i = 0; i < calibrationStats.size(); i++){
            timePerc = (calibrationStats.at(i).duration /
                        totalDurationMs) * 100.0;

            _calibrationFile << calibrationStats.at(i).numSteps << "\t";
            _calibrationFile << calibrationStats.at(i).duration << "\t";
            _calibrationFile << timePerc << "\t";
            _calibrationFile << std::endl;
        }
    }
};

/*!
 * \class adp_ff_farm
 * \brief This class wraps a farm to let it reconfigurable.
 *
 * This class wraps a farm to let it reconfigurable.
 */
template<typename lb_t=ff_loadbalancer, typename gt_t=ff_gatherer>
class adp_ff_farm: public ff_farm<lb_t, gt_t>{
    friend class AdaptivityManagerFarm;
private:
    std::vector<adp_ff_node*> _adaptiveWorkers;
    adp_ff_node* _adaptiveEmitter;
    adp_ff_node* _adaptiveCollector;
    bool _firstRun;
    AdaptivityParameters _adaptivityParameters;
    AdaptivityManagerFarm* _adaptivityManager;

    void construct();

    void construct(AdaptivityParameters adaptivityParameters);

    std::vector<adp_ff_node*> getAdaptiveWorkers() const;

    adp_ff_node* getAdaptiveEmitter() const;

    adp_ff_node* getAdaptiveCollector() const;

    void waitInternal();
public:

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     */
    adp_ff_farm(std::vector<ff_node*>& w, ff_node* const emitter = NULL,
                ff_node* const collector = NULL, bool inputCh = false);

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     */
    adp_ff_farm(bool inputCh = false,
                int inBufferEntries = ff_farm<lb_t, gt_t>::DEF_IN_BUFF_ENTRIES,
                int outBufferEntries = ff_farm<lb_t, gt_t>::DEF_OUT_BUFF_ENTRIES,
                bool workerCleanup = false,
                int maxNumWorkers = DEF_MAX_NUM_WORKERS,
                bool fixedSize = false);

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     * @param adaptivityParameters Parameters that will be used by the farm
     *                             to take reconfiguration decisions.
     */
    adp_ff_farm(AdaptivityParameters adaptivityParameters,
                std::vector<ff_node*>& w,
                ff_node* const emitter = NULL,
                ff_node* const collector = NULL,
                bool inputCh = false);

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     * @param adaptivityParameters Parameters that will be used by the farm
     *                             to take reconfiguration decisions.
     */
    explicit
    adp_ff_farm(AdaptivityParameters adaptivityParameters,
                bool inputCh = false,
                int inBufferEntries = ff_farm<lb_t, gt_t>::DEF_IN_BUFF_ENTRIES,
                int outBufferEntries = ff_farm<lb_t, gt_t>::DEF_OUT_BUFF_ENTRIES,
                bool workerCleanup = false,
                int maxNumWorkers = DEF_MAX_NUM_WORKERS,
                bool fixedSize = false);

    /**
     * Destroyes this adaptive farm.
     */
    ~adp_ff_farm();

    void setAdaptivityParameters(AdaptivityParameters adaptivityParameters);

    void firstRunBefore();

    void firstRunAfter();

    /**
     * Runs this farm.
     */
    int run(bool skip_init=false);

    /**
     * Waits this farm for completion.
     */
    int wait();
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
    cpufreq::Frequency frequency;

    FarmConfiguration():numWorkers(0), frequency(0){;}
    FarmConfiguration(uint numWorkers, cpufreq::Frequency frequency = 0):
                     numWorkers(numWorkers), frequency(frequency){;}
}FarmConfiguration;

std::ostream& operator<<(std::ostream& os, const FarmConfiguration& obj){
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
    energy::JoulesCpu watts; ///< Watts consumed by all the CPUs.
    double bandwidth; ///< Bandwidth of the entire farm.
    double utilization; ///< Utilization of the entire farm.

    MonitoredSample():bandwidth(0), utilization(0){;}

    void swap(MonitoredSample& x){
        using std::swap;

        swap(watts, x.watts);
        swap(bandwidth, x.bandwidth);
        swap(utilization, x.utilization);
    }

    MonitoredSample& operator=(MonitoredSample rhs){
        swap(rhs);
        return *this;
    }

    MonitoredSample& operator+=(const MonitoredSample& rhs){
        watts += rhs.watts;
        bandwidth += rhs.bandwidth;
        utilization += rhs.utilization;
        return *this;
    }

    MonitoredSample& operator-=(const MonitoredSample& rhs){
        watts -= rhs.watts;
        bandwidth -= rhs.bandwidth;
        utilization -= rhs.utilization;
        return *this;
    }

    MonitoredSample& operator*=(const MonitoredSample& rhs){
        watts *= rhs.watts;
        bandwidth *= rhs.bandwidth;
        utilization *= rhs.utilization;
        return *this;
    }

    MonitoredSample& operator/=(const MonitoredSample& rhs){
        watts /= rhs.watts;
        bandwidth /= rhs.bandwidth;
        utilization /= rhs.utilization;
        return *this;
    }

    MonitoredSample operator/=(double x){
        watts /= x;
        bandwidth /= x;
        utilization /= x;
        return *this;
    }

    MonitoredSample operator*=(double x){
        watts *= x;
        bandwidth *= x;
        utilization *= x;
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

std::ostream& operator<<(std::ostream& os, const MonitoredSample& obj){
    os << "[";
    os << "Watts: " << obj.watts << " ";
    os << "Bandwidth: " << obj.bandwidth << " ";
    os << "Utilization: " << obj.utilization << " ";
    os << "]";
    return os;
}

std::ofstream& operator<<(std::ofstream& os, const MonitoredSample& obj){
    os << obj.watts.cores << "\t";
    os << obj.bandwidth << "\t";
    os << obj.utilization << "\t";
    return os;
}

inline energy::JoulesCpu squareRoot(const energy::JoulesCpu& x){
    energy::JoulesCpu r;
    r.cores = x.cores?std::sqrt(x.cores):0;
    r.cpu = x.cpu?std::sqrt(x.cpu):0;
    r.graphic = x.graphic?std::sqrt(x.graphic):0;
    r.dram = x.dram?std::sqrt(x.dram):0;
    return r;
}

inline MonitoredSample squareRoot(const MonitoredSample& x){
    MonitoredSample r;
    r.watts = squareRoot(x.watts);
    r.bandwidth = x.bandwidth?std::sqrt(x.bandwidth):0;
    r.utilization = x.utilization?std::sqrt(x.utilization):0;
    return r;
}


inline void regularize(MonitoredSample& x){
    if(x.watts.cpu < 0){x.watts.cpu = 0;}
    if(x.watts.cores < 0){x.watts.cores = 0;}
    if(x.watts.graphic < 0){x.watts.graphic = 0;}
    if(x.watts.dram < 0){x.watts.dram = 0;}
    if(x.bandwidth < 0){x.bandwidth = 0;}
    if(x.utilization < 0){x.utilization = 0;}
}

/*!
 * \internal
 * \class AdaptivityManagerFarm
 * \brief This class manages the adaptivity in farm based computations.
 *
 * This class manages the adaptivity in farm based computations.
 */
class AdaptivityManagerFarm: public utils::Thread{
    friend class PredictorSimple;
    friend class PredictorLinearRegression;
    friend class RegressionData;
    friend class RegressionDataServiceTime;
    friend class RegressionDataPower;
    friend class Calibrator;
    friend class CalibratorLowDiscrepancy;
private:
    utils::Monitor _monitor; ///< Used to let the manager stop safely.
    adp_ff_farm<>* _farm; ///< The managed farm.
    AdaptivityParameters _p; ///< The parameters used to take management decisions.
    uint _startTimeMs; ///< Starting time of the manager.
    cpufreq::CpuFreq* _cpufreq; ///< The cpufreq module.
    energy::Energy* _energy; ///< The energy module.
    task::TasksManager* _task; ///< The task module.
    topology::Topology* _topology; ///< The topology module.
    uint _numCpus; ///< Number of CPUs
    uint _numPhysicalCores; ///< Number of physical cores.
    uint _numPhysicalCoresPerCpu; ///< Number of physical cores per CPU.
    uint _numVirtualCoresPerPhysicalCore; ///< Number of virtual cores per physical core.
    adp_ff_node* _emitter; ///< The emitter (if present).
    adp_ff_node* _collector; ///< The collector (if present).
    std::vector<adp_ff_node*> _activeWorkers; ///< The currently running workers.
    std::vector<adp_ff_node*> _inactiveWorkers; ///< Workers that can run but are not currently running.
    size_t _maxNumWorkers; ///< The maximum number of workers that can be activated by the manager.
    bool _emitterSensitivitySatisfied; ///< If true, the user requested sensitivity for emitter and the
                                       ///< request has been satisfied.
    bool _collectorSensitivitySatisfied; ///< If true, the user requested sensitivity for collector and the
                                         ///< request has been satisfied.
    FarmConfiguration _currentConfiguration; ///< The current configuration of the farm.
    std::vector<topology::CpuId> _usedCpus; ///< CPUs currently used by farm nodes.
    std::vector<topology::CpuId> _unusedCpus; ///< CPUs not used by farm nodes.
    const std::vector<topology::VirtualCore*> _availableVirtualCores; ///< The available virtual cores, sorted according to
                                                                      ///< the mapping strategy.
    std::vector<topology::VirtualCore*> _activeWorkersVirtualCores; ///< The virtual cores where the active workers
                                                                    ///< are running.
    std::vector<topology::VirtualCore*> _inactiveWorkersVirtualCores; ///< The virtual cores where the inactive workers
                                                                      ///< are running.
    std::vector<topology::VirtualCore*> _unusedVirtualCores; ///< Virtual cores not used by the farm nodes.
    topology::VirtualCore* _emitterVirtualCore; ///< The virtual core where the emitter (if present) is running.
    topology::VirtualCore* _collectorVirtualCore; ///< The virtual core where the collector (if present) is running.
    std::vector<cpufreq::Domain*> _scalableDomains; ///< The domains on which frequency scaling is applied.
    cpufreq::VoltageTable _voltageTable; ///< The voltage table.
    std::vector<cpufreq::Frequency> _availableFrequencies; ///< The available frequencies on this machine.
    MovingAverage<MonitoredSample>* _monitoredSamples; ///< Monitored samples;
    double _totalTasks; ///< The number of tasks processed since the last reconfiguration.
    double _averageBandwidth; ///< The average tasks per second processed during last time window.
    double _averageUtilization; ///< The average utilization during last time window.
    energy::JoulesCpu _averageWatts; ///< Average Watts during last time window.
    uint64_t _remainingTasks; ///< When contract is CONTRACT_COMPLETION_TIME, represent the number of tasks that
                             ///< still needs to be processed by the application.
    time_t _deadline; ///< When contract is CONTRACT_COMPLETION_TIME, represent the deadline of the application.
    double _lastStoredSampleMs; ///< Milliseconds timestamp of the last store of a sample.

    Calibrator* _calibrator; ///< The calibrator of the predictors.
    Predictor* _primaryPredictor; ///< The predictor of the primary value.
    Predictor* _secondaryPredictor; ///< The predictor of the secondary value.
    double _primaryPrediction; ///< The prediction done for the primary value for the chosen configuration.
    double _secondaryPrediction; ///< The prediction done for the secondary value for the chosen configuration.

#ifdef DEBUG_FARM
    std::ofstream samplesFile;
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
     * If possible, finds a set of physical cores belonging to domains different from
     * those of virtual cores in 'virtualCores' vector.
     * @param virtualCores A vector of virtual cores.
     * @return A set of physical cores that can always run at the highest frequency.
     */
    std::vector<topology::PhysicalCore*> getSeparatedDomainPhysicalCores(const std::vector<topology::VirtualCore*>& virtualCores) const{
        std::vector<cpufreq::Domain*> allDomains = _cpufreq->getDomains();
        std::vector<cpufreq::Domain*> hypotheticWorkersDomains = _cpufreq->getDomains(virtualCores);
        std::vector<topology::PhysicalCore*> physicalCoresInUnusedDomains;
        if(allDomains.size() > hypotheticWorkersDomains.size()){
           for(size_t i = 0; i < allDomains.size(); i++){
               cpufreq::Domain* currentDomain = allDomains.at(i);
               if(!utils::contains(hypotheticWorkersDomains, currentDomain)){
                   utils::insertToEnd(_topology->virtualToPhysical(currentDomain->getVirtualCores()),
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
    void setDomainToHighestFrequency(const cpufreq::Domain* domain){
        if(!domain->setGovernor(mammut::cpufreq::GOVERNOR_PERFORMANCE)){
            if(!domain->setGovernor(mammut::cpufreq::GOVERNOR_USERSPACE) ||
               !domain->setHighestFrequencyUserspace()){
                throw std::runtime_error("AdaptivityManagerFarm: Fatal error while setting highest frequency for "
                                         "sensitive emitter/collector. Try to run it without sensitivity parameters.");
            }
        }
    }

    /**
     * Computes the available virtual cores, sorting them according to the specified
     * mapping strategy.
     * @return The available virtual cores, sorted according to the specified mapping
     *         strategy.
     */
    std::vector<topology::VirtualCore*> getAvailableVirtualCores(){
        std::vector<topology::VirtualCore*> r;
        if(_p.strategyMapping == STRATEGY_MAPPING_AUTO){
            _p.strategyMapping = STRATEGY_MAPPING_LINEAR;
        }

        switch(_p.strategyMapping){
            case STRATEGY_MAPPING_LINEAR:{
               /*
                * Generates a vector of virtual cores to be used for linear mapping.node
                * It contains first one virtual core per physical core (virtual cores
                * on the same CPU are consecutive).
                * Then, the other groups of virtual cores follow.
                */
                std::vector<topology::Cpu*> cpus = _topology->getCpus();

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
                        std::vector<topology::PhysicalCore*> physicalCores = cpus.at(i)->getPhysicalCores();
                        for(size_t j = 0; j < physicalCores.size(); j++){
                            r.push_back(physicalCores.at(j)->getVirtualCores().at(virtualUsed));
                        }
                    }
                    ++virtualUsed;
                }
            }break;
            case STRATEGY_MAPPING_CACHE_EFFICIENT:{
                throw std::runtime_error("Not yet supported.");
            }
            default:
                break;
        }
        return r;
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

        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO &&
          ((_p.mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE && !_emitterSensitivitySatisfied) ||
           (_p.mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE && !_collectorSensitivitySatisfied))){
            size_t scalableVirtualCoresNum = _activeWorkers.size() +
                                             (_emitter && _p.mappingEmitter != SERVICE_NODE_MAPPING_PERFORMANCE) +
                                             (_collector && _p.mappingEmitter != SERVICE_NODE_MAPPING_PERFORMANCE);
            /** When sensitive is specified, we always choose the WEC mapping. **/
            std::vector<topology::VirtualCore*> scalableVirtualCores(_availableVirtualCores.begin(), (scalableVirtualCoresNum < _availableVirtualCores.size())?
                                                                                                   _availableVirtualCores.begin() + scalableVirtualCoresNum:
                                                                                                   _availableVirtualCores.end());
            std::vector<topology::PhysicalCore*> performancePhysicalCores;
            performancePhysicalCores = getSeparatedDomainPhysicalCores(scalableVirtualCores);
            if(performancePhysicalCores.size()){
                size_t index = 0;

                if(_p.mappingEmitter == SERVICE_NODE_MAPPING_PERFORMANCE){
                    topology::VirtualCore* vc = performancePhysicalCores.at(index)->getVirtualCore();
                    setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                    _emitterVirtualCore = vc;
                    _emitterSensitivitySatisfied = true;
                    index = (index + 1) % performancePhysicalCores.size();
                }

                if(_p.mappingCollector == SERVICE_NODE_MAPPING_PERFORMANCE){
                    topology::VirtualCore* vc = performancePhysicalCores.at(index)->getVirtualCore();
                    setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                    _collectorVirtualCore = vc;
                    _collectorSensitivitySatisfied = true;
                    index = (index + 1) % performancePhysicalCores.size();
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
        nextIndex = (nextIndex + _activeWorkers.size()) % _availableVirtualCores.size();

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
        //      adaptive node e quindi la getThreadHandler non c'e' (fallisce lo static_cast prima)):(
        if(_emitter){
            if(!_emitterVirtualCore){
                _emitterVirtualCore = _availableVirtualCores.at(emitterIndex);
            }
            if(!_emitterVirtualCore->isHotPlugged()){
                _emitterVirtualCore->hotPlug();
            }
            _emitter->getThreadHandler()->move(_emitterVirtualCore);
        }

        for(size_t i = 0; i < _activeWorkers.size(); i++){
            topology::VirtualCore* vc = _availableVirtualCores.at((firstWorkerIndex + i) % _availableVirtualCores.size());
            _activeWorkersVirtualCores.push_back(vc);
            if(!vc->isHotPlugged()){
                vc->hotPlug();
            }
            _activeWorkers.at(i)->getThreadHandler()->move(vc);
        }

        if(_collector){
            if(!_collectorVirtualCore){
                _collectorVirtualCore = _availableVirtualCores.at(collectorIndex);
            }
            if(!_collectorVirtualCore->isHotPlugged()){
                _collectorVirtualCore->hotPlug();
            }
            _collector->getThreadHandler()->move(_collectorVirtualCore);
        }
    }

    /**
     * Apply a specific strategy for a specified set of virtual cores.
     * @param strategyUnused The unused strategy.
     * @param unusedVirtualCores The virtual cores.
     */
    void applyUnusedVirtualCoresStrategy(StrategyUnusedVirtualCores strategyUnused,
                                         const std::vector<topology::VirtualCore*>& unusedVirtualCores){
        switch(strategyUnused){
            case STRATEGY_UNUSED_VC_OFF:{
                for(size_t i = 0; i < unusedVirtualCores.size(); i++){
                    topology::VirtualCore* vc = unusedVirtualCores.at(i);
                    if(vc->isHotPluggable() && vc->isHotPlugged()){
                        vc->hotUnplug();
                    }
                }
            }break;
            case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                std::vector<cpufreq::Domain*> unusedDomains = _cpufreq->getDomainsComplete(unusedVirtualCores);
                for(size_t i = 0; i < unusedDomains.size(); i++){
                    cpufreq::Domain* domain = unusedDomains.at(i);
                    if(!domain->setGovernor(cpufreq::GOVERNOR_POWERSAVE)){
                        if(!domain->setGovernor(cpufreq::GOVERNOR_USERSPACE) ||
                           !domain->setLowestFrequencyUserspace()){
                            throw std::runtime_error("AdaptivityManagerFarm: Impossible to set lowest frequency "
                                                     "for unused virtual cores.");
                        }
                    }
                }
            }break;
            default:{
                return;
            }
        }
    }

    /**
     * Apply the strategies for inactive and unused virtual cores.
     */
    void applyUnusedVirtualCoresStrategy(){
        /**
         * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the virtual cores
         * on a domain, we can also lower its frequency to the minimum.
         */
        std::vector<topology::VirtualCore*> virtualCores;
        if(_p.strategyInactiveVirtualCores != STRATEGY_UNUSED_VC_NONE){
            utils::insertToEnd(_inactiveWorkersVirtualCores, virtualCores);
        }
        if(_p.strategyUnusedVirtualCores != STRATEGY_UNUSED_VC_NONE){
            utils::insertToEnd(_unusedVirtualCores, virtualCores);
        }
        applyUnusedVirtualCoresStrategy(STRATEGY_UNUSED_VC_LOWEST_FREQUENCY, virtualCores);


        virtualCores.clear();
        if(_p.strategyInactiveVirtualCores == STRATEGY_UNUSED_VC_OFF ||
           _p.strategyInactiveVirtualCores == STRATEGY_UNUSED_VC_AUTO){
            utils::insertToEnd(_inactiveWorkersVirtualCores, virtualCores);
        }
        if(_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_OFF ||
           _p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_AUTO){
            utils::insertToEnd(_unusedVirtualCores, virtualCores);
        }
        applyUnusedVirtualCoresStrategy(STRATEGY_UNUSED_VC_AUTO, virtualCores);
    }

    /**
     * Updates the scalable domains vector.
     */
    void updateScalableDomains(){
        std::vector<topology::VirtualCore*> frequencyScalableVirtualCores = _activeWorkersVirtualCores;
        /**
         * Node sensitivity may be not satisfied both because it was not requested, or because
         * it was requested but it was not possible to satisfy it. In both cases, we need to scale
         * the virtual core of the node as we scale the others.
         */
        if(_emitter && !_emitterSensitivitySatisfied){
            frequencyScalableVirtualCores.push_back(_emitterVirtualCore);
        }
        if(_collector && !_collectorSensitivitySatisfied){
            frequencyScalableVirtualCores.push_back(_collectorVirtualCore);
        }

        _scalableDomains = _cpufreq->getDomains(frequencyScalableVirtualCores);
    }

    /**
     * Set a specific P-state for the virtual cores used by
     * the current active workers, emitter (if not sensitive) and
     * collector (if not sensitive).
     * @param frequency The frequency to be set.
     */
    void updatePstate(cpufreq::Frequency frequency){
        updateScalableDomains();
        cpufreq::Domain* currentDomain;
        for(size_t i = 0; i < _scalableDomains.size(); i++){
            currentDomain = _scalableDomains.at(i);
            if(!currentDomain->setGovernor(_p.frequencyGovernor)){
                throw std::runtime_error("AdaptivityManagerFarm: Impossible to "
                                         "set the specified governor.");
            }
            if(_p.frequencyGovernor != cpufreq::GOVERNOR_USERSPACE){
                if(!currentDomain->setGovernorBounds(_p.frequencyLowerBound,
                                                     _p.frequencyUpperBound)){
                    throw std::runtime_error("AdaptivityManagerFarm: Impossible "
                                             "to set the specified governor's "
                                             "bounds.");
                }
            }else if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                if(!currentDomain->setFrequencyUserspace(frequency)){
                    throw std::runtime_error("AdaptivityManagerFarm: Impossible "
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
            topology::VirtualCore* vc = _availableVirtualCores.at(i);
            if(vc != _emitterVirtualCore && vc != _collectorVirtualCore &&
               !utils::contains(_activeWorkersVirtualCores, vc) &&
               !utils::contains(_inactiveWorkersVirtualCores, vc)){
                _unusedVirtualCores.push_back(vc);
            }
        }
        updateUsedCpus();
        applyUnusedVirtualCoresStrategy();

        // Insert dummy constant frequency
        _availableFrequencies.push_back(1.0);

        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO){
            if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                // We suppose that all the domains have the same
                // available frequencies.
                _availableFrequencies = _cpufreq->getDomains().at(0)->getAvailableFrequencies();

                // Remove turbo boost frequency.
                if(!_p.turboBoost){
                    if(utils::intToString(_availableFrequencies.back()).at(3) == '1'){
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
     * Updates the monitored values.
     */
    void updateMonitoredValues(){
        MonitoredSample avg = _monitoredSamples->average();
        _averageBandwidth = avg.bandwidth;
        _averageUtilization = avg.utilization;
        _averageWatts = avg.watts;
    }

    /**
     * Returns the primary value according to the required contract.
     * @return The primary value according to the required contract.
     */
    double getPrimaryValue() const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return _averageUtilization;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return _averageBandwidth;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return _averageWatts.cores;
            }break;
            default:{
                return 0;
            }break;
        }
    }

    /**
     * Returns the secondary value according to the required contract.
     * @return The secondary value according to the required contract.
     */
    double getSecondaryValue() const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return _averageWatts.cores;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return _averageBandwidth;
            }break;
            default:{
                return 0;
            }break;
        }
    }

    /**
     * Checks if the contract is violated.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated() const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return _averageUtilization < _p.underloadThresholdFarm ||
                       _averageUtilization > _p.overloadThresholdFarm;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                double tolerance = (_p.requiredBandwidth *
                                    _p.maxPredictionError) / 100.0;
                return _averageBandwidth < _p.requiredBandwidth - tolerance;
            }break;
            case CONTRACT_POWER_BUDGET:{
                double tolerance = (_p.powerBudget *
                                    _p.maxPredictionError) / 100.0;
                return _averageWatts.cores > _p.powerBudget + tolerance;
            }break;
            default:{
                return false;
            }break;
        }
        return false;
    }

    /**
     * Checks if a specific primary value is feasible according to the
     * specified contract.
     * @return true if the solution is feasible, false otherwise.
     */
    bool isFeasiblePrimaryValue(double primaryValue) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return primaryValue > _p.underloadThresholdFarm &&
                        primaryValue < _p.overloadThresholdFarm;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return primaryValue > _p.requiredBandwidth;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return primaryValue < _p.powerBudget;
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
        cpufreq::VoltageTableKey key(configuration.numWorkers,
                                     configuration.frequency);
        cpufreq::VoltageTableIterator it = _voltageTable.find(key);
        if(it != _voltageTable.end()){
            return it->second;
        }else{
            throw std::runtime_error("Frequency and/or number of virtual cores "
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
                    return std::abs(distanceX) < std::abs(distanceY);
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
        FarmConfiguration r;
        double currentPrimaryPrediction = 0;
        double currentSecondaryPrediction = 0;

        double primaryPrediction = 0;
        double secondaryPrediction = 0;
        double primaryPredictionSub = 0;
        double secondaryPredictionSub = 0;

        double bestSuboptimalValue = getPrimaryValue();
        FarmConfiguration bestSuboptimalConfiguration = _currentConfiguration;
        bool feasibleSolutionFound = false;

        double bestSecondaryValue = 0;
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                // We have to minimize the power/energy.
                bestSecondaryValue = std::numeric_limits<double>::max();
            }break;
            case CONTRACT_POWER_BUDGET:{
                // We have to maximize the bandwidth.
                bestSecondaryValue = std::numeric_limits<double>::min();
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
                FarmConfiguration examinedConfiguration(i, _availableFrequencies.at(j));
                currentPrimaryPrediction = _primaryPredictor->predict(examinedConfiguration);
                switch(_p.contractType){
                    case CONTRACT_PERF_COMPLETION_TIME:{
                        remainingTime = (double) _remainingTasks / currentPrimaryPrediction;
                    }break;
                    case CONTRACT_PERF_UTILIZATION:{
                        currentPrimaryPrediction = (_averageBandwidth / currentPrimaryPrediction) * _averageUtilization;
                    }break;
                    default:{
                        ;
                    }
                }

                if(isFeasiblePrimaryValue(currentPrimaryPrediction)){
                    currentSecondaryPrediction = _secondaryPredictor->predict(examinedConfiguration);
                    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                        currentSecondaryPrediction *= remainingTime;
                    }
                    if(isBestSecondaryValue(currentSecondaryPrediction,
                                            bestSecondaryValue)){
                        bestSecondaryValue = currentSecondaryPrediction;
                        r = examinedConfiguration;
                        feasibleSolutionFound = true;
                        primaryPrediction = currentPrimaryPrediction;
                        secondaryPrediction = currentSecondaryPrediction;
                    }
                }else if(!feasibleSolutionFound &&
                         isBestSuboptimalValue(currentPrimaryPrediction,
                                               bestSuboptimalValue)){
                    bestSuboptimalValue = currentPrimaryPrediction;
                    bestSuboptimalConfiguration = examinedConfiguration;
                    primaryPredictionSub = currentPrimaryPrediction;
                    secondaryPredictionSub = currentSecondaryPrediction;
                }
            }
        }

        if(feasibleSolutionFound){
            _primaryPrediction = primaryPrediction;
            _secondaryPrediction = secondaryPrediction;
            return r;
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
            topology::CpuId cpuId = _activeWorkersVirtualCores.at(i)->getCpuId();
            if(!utils::contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }
        if(_emitterVirtualCore){
            topology::CpuId cpuId = _emitterVirtualCore->getCpuId();
            if(!utils::contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }
        if(_collectorVirtualCore){
            topology::CpuId cpuId = _collectorVirtualCore->getCpuId();
            if(!utils::contains(_usedCpus, cpuId)){
                _usedCpus.push_back(cpuId);
            }
        }

        std::vector<topology::Cpu*> cpus = _topology->getCpus();
        for(size_t i = 0; i < cpus.size(); i++){
            if(!utils::contains(_usedCpus, cpus.at(i)->getCpuId())){
                _unusedCpus.push_back(cpus.at(i)->getCpuId());
            }
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
            throw std::runtime_error("AdaptivityManagerFarm: fatal error, trying to activate zero "
                                     "workers.");
        }else if(configuration.numWorkers > _maxNumWorkers){
            throw std::runtime_error("AdaptivityManagerFarm: fatal error, trying to activate more "
                                     "workers than the maximum allowed.");
        }

        /****************** Refine the model ******************/
        if(refine){
            _primaryPredictor->refine();
            _secondaryPredictor->refine();
        }

        /****************** Workers change started ******************/
        if(_currentConfiguration.numWorkers != configuration.numWorkers){
            std::vector<cpufreq::RollbackPoint> rollbackPoints;
            if(_p.fastReconfiguration){
                for(size_t i = 0; i < _scalableDomains.size(); i++){
                    rollbackPoints.push_back(_scalableDomains.at(i)->getRollbackPoint());
                    setDomainToHighestFrequency(_scalableDomains.at(i));
                }
            }

            /** Stops farm. **/
            DEBUG("Asking the farm for NULL production");
            _emitter->produceNull();
            DEBUG("Wait for freezing");
            _farm->wait_freezing();
            DEBUG("Farm freezed");
            /**
             * When workers stops, they update their last sample.
             * Accordingly, we do not need to ask them but we can just
             * retrieve the samples.
             */
            storeNewSample(false);

            if(_currentConfiguration.numWorkers > configuration.numWorkers){
                /** Move workers from active to inactive. **/
                uint workersNumDiff = _currentConfiguration.numWorkers - configuration.numWorkers;
                utils::moveEndToFront(_activeWorkers, _inactiveWorkers, workersNumDiff);
                utils::moveEndToFront(_activeWorkersVirtualCores, _inactiveWorkersVirtualCores, workersNumDiff);
            }else{
                /** Move workers from inactive to active. **/
                /**
                 * We need to map them again because if virtual cores were shutdown, the
                 * threads that were running on them have been moved on different virtual
                 * cores. Accordingly, when started again they would not be run on the
                 * correct virtual cores.
                 */
                uint workersNumDiff = configuration.numWorkers - _currentConfiguration.numWorkers;
                for(size_t i = 0; i < workersNumDiff; i++){
                    topology::VirtualCore* vc = _inactiveWorkersVirtualCores.at(i);
                    if(!vc->isHotPlugged()){
                        vc->hotPlug();
                    }
                    _inactiveWorkers.at(i)->getThreadHandler()->move(vc);
                }
                utils::moveFrontToEnd(_inactiveWorkers, _activeWorkers, workersNumDiff);
                utils::moveFrontToEnd(_inactiveWorkersVirtualCores, _activeWorkersVirtualCores, workersNumDiff);
            }

            updateUsedCpus();

            /** Notify the nodes that a reconfiguration is happening. **/
            _emitter->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            for(uint i = 0; i < configuration.numWorkers; i++){
                _activeWorkers.at(i)->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            }
            if(_collector){
                _collector->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            }


            /** Start the farm again. **/
            DEBUG("Asking the farm to start again");
            _farm->run_then_freeze(configuration.numWorkers);
            DEBUG("Farm started");
            //TODO: Se la farm non è stata avviata con la run_then_freeze questo potrebbe essere un problema.
            if(_p.fastReconfiguration){
                _cpufreq->rollback(rollbackPoints);
            }
        }
        /****************** Workers change terminated ******************/
        applyUnusedVirtualCoresStrategy();
        /****************** P-state change started ******************/
        //TODO: Maybe sensitivity could not be satisfied with the maximum number of workers but could be satisfied now.
        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO){
            updatePstate(configuration.frequency);
        }
        /****************** P-state change terminated ******************/
        _currentConfiguration = configuration;

        /****************** Clean state ******************/
        _lastStoredSampleMs = utils::getMillisecondsTime();
        _monitoredSamples->reset();
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
                                 _monitoredSamples->getLastSample().bandwidth,
                                 _averageBandwidth,
                                 _monitoredSamples->coefficientVariation().bandwidth,
                                 _averageUtilization,
                                 _averageWatts);
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
        sample = WorkerSample();
        for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
            WorkerSample tmp;
            if(!_activeWorkers.at(i)->getSampleResponse(tmp)){
                return false;
            }
            sample += tmp;
        }
        sample.loadPercentage /= _currentConfiguration.numWorkers;
        return true;
    }


    /**
     * Store a new sample.
     * @param askWorkers If false, the request to the worker is not sent.
     * @return false if the farm is not running anymore, true otherwise.
     **/
    bool storeNewSample(bool askWorkers = true){
        if(askWorkers){
            askForWorkersSamples();
        }

        MonitoredSample sample;
        WorkerSample ws;
        energy::JoulesCpu joules;
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
            energy::CounterCpu* currentCounter = _energy->getCounterCpu(_usedCpus.at(i));
            joules += currentCounter->getJoules();
        }
        for(size_t i = 0; i < _unusedCpus.size(); i++){
            energy::CounterCpu* currentCounter = _energy->getCounterCpu(_unusedCpus.at(i));
            joules += currentCounter->getJoules();
        }

        double now = utils::getMillisecondsTime();
        double durationSecs = (now - _lastStoredSampleMs) / 1000.0;
        _lastStoredSampleMs = now;

        sample.watts = joules / durationSecs;
        sample.utilization = ws.loadPercentage;
        sample.bandwidth = (double) ws.tasksCount / durationSecs;

        _energy->resetCountersCpu();
        _monitoredSamples->add(sample);

        DEBUGB(samplesFile << *_monitoredSamples << "\n");
        return true;
    }


    void initCalibrators(){
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
                _primaryPredictor = new PredictorLinearRegression(primary, *this);
                _secondaryPredictor = new PredictorLinearRegression(secondary, *this);
                initCalibrators();
            }break;
            default:{
                ;
            }break;
        }
    }

public:
    /**
     * Creates a farm adaptivity manager.
     * ATTENTION: Needs to be created when the farm is ready (i.e. in run* methods).
     * @param farm The farm to be managed.
     * @param adaptivityParameters The parameters to be used for adaptivity decisions.
     */
    AdaptivityManagerFarm(adp_ff_farm<>* farm, AdaptivityParameters adaptivityParameters):
        _farm(farm),
        _p(adaptivityParameters),
        _startTimeMs(0),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _energy(_p.mammut.getInstanceEnergy()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _numCpus(_topology->getCpus().size()),
        _numPhysicalCores(_topology->getPhysicalCores().size()),
        _numPhysicalCoresPerCpu(_topology->getCpu(0)->getPhysicalCores().size()),
        _numVirtualCoresPerPhysicalCore(_topology->getPhysicalCore(0)->getVirtualCores().size()),
        _emitter(_farm->getAdaptiveEmitter()),
        _collector(_farm->getAdaptiveCollector()),
        _activeWorkers(_farm->getAdaptiveWorkers()),
        _maxNumWorkers(_activeWorkers.size()),
        _emitterSensitivitySatisfied(false),
        _collectorSensitivitySatisfied(false),
        _currentConfiguration(_activeWorkers.size(), 0),
        _availableVirtualCores(getAvailableVirtualCores()),
        _emitterVirtualCore(NULL),
        _collectorVirtualCore(NULL){
        /** If voltage table file is specified, then load the table. **/
        if(_p.voltageTableFile.compare("")){
            cpufreq::loadVoltageTable(_voltageTable, _p.voltageTableFile);
        }

        _monitoredSamples = NULL;
        switch(_p.strategySmoothing){
        case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
            _monitoredSamples = new MovingAverageSimple<MonitoredSample>
                                    (_p.numSamples);
        }break;
        case STRATEGY_SMOOTHING_EXPONENTIAL:{
            _monitoredSamples = new MovingAverageExponential<MonitoredSample>
                                    (_p.alphaExpAverage);
        }break;
        }

        DEBUGB(samplesFile.open("samples.csv"));
    }

    /**
     * Destroyes this adaptivity manager.
     */
    ~AdaptivityManagerFarm(){
        delete _monitoredSamples;
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
        /**
         * Wait for all the nodes to be running.
         */
        for(size_t i = 0; i < _activeWorkers.size(); i++){
            _activeWorkers.at(i)->waitThreadCreation();
        }
        if(_emitter){
            _emitter->waitThreadCreation();
        }
        if(_collector){
            _collector->waitThreadCreation();
        }

        _startTimeMs = utils::getMillisecondsTime();

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
            _deadline = time(NULL) + _p.requiredCompletionTime;
        }

        initPredictors();

        double microsecsSleep = 0;
        if(_p.contractType == CONTRACT_NONE){
            //_monitor.wait();
            _farm->waitInternal();
            storeNewSample();
            updateMonitoredValues();
            observe();
        }else{
            /* Force the first calibration point. **/
            if(_calibrator){
                changeConfiguration(_calibrator->getNextConfiguration(), false);
            }

            double startSample = utils::getMillisecondsTime();
            while(!mustStop()){
                double overheadMs = utils::getMillisecondsTime() - startSample;
                microsecsSleep = ((double)_p.samplingInterval - overheadMs)*
                                  (double)MAMMUT_MICROSECS_IN_MILLISEC;
                if(microsecsSleep < 0){
                    microsecsSleep = 0;
                }
                usleep(microsecsSleep);
                startSample = utils::getMillisecondsTime();

                if(!storeNewSample()){
                    goto controlLoopEnd;
                }

                if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                    time_t now = time(NULL);
                    if(now >= _deadline){
                        _p.requiredBandwidth = std::numeric_limits<double>::max();
                    }else{
                        _p.requiredBandwidth = _remainingTasks / (_deadline - now);
                    }
                }

                updateMonitoredValues();
                observe();

                if(_monitoredSamples->size() >= _p.numSamples){
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
                        startSample = utils::getMillisecondsTime();
                    }
                }
            }
        }
    controlLoopEnd:

        if(_calibrator && _p.observer){
            _p.observer->calibrationStats(_calibrator->getCalibrationsStats(),
                                          utils::getMillisecondsTime() -
                                          _startTimeMs);
        }
    }

    /**
     * Stops this manager.
     */
    void stop(){
        _monitor.notifyOne();
    }

    /**
     * Return true if the manager must stop, false otherwise.
     * @retur true if the manager must stop, false otherwise.
     */
    bool mustStop(){
        return _monitor.predicate();
    }
};

template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::construct(){
    _adaptiveEmitter = NULL;
    _adaptiveCollector = NULL;
    _firstRun = true;
    _adaptivityManager = NULL;
}

template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::construct(AdaptivityParameters adaptivityParameters){
    construct();
    setAdaptivityParameters(adaptivityParameters);
}

template <typename lb_t, typename gt_t>
std::vector<adp_ff_node*> adp_ff_farm<lb_t, gt_t>::getAdaptiveWorkers() const{
    return _adaptiveWorkers;
}

template <typename lb_t, typename gt_t>
adp_ff_node* adp_ff_farm<lb_t, gt_t>::getAdaptiveEmitter() const{
    return _adaptiveEmitter;
}

template <typename lb_t, typename gt_t>
adp_ff_node* adp_ff_farm<lb_t, gt_t>::getAdaptiveCollector() const{
    return _adaptiveCollector;
}

template <typename lb_t, typename gt_t>
adp_ff_farm<lb_t, gt_t>::adp_ff_farm(std::vector<ff_node*>& w, ff_node* const emitter,
                                     ff_node* const collector, bool inputCh):
    ff_farm<lb_t, gt_t>::ff_farm(w, emitter, collector, inputCh){
    construct();
}

template <typename lb_t, typename gt_t>
adp_ff_farm<lb_t, gt_t>::adp_ff_farm(bool inputCh,
                      int inBufferEntries,
                      int outBufferEntries,
                      bool workerCleanup,
                      int maxNumWorkers,
                      bool fixedSize):
    ff_farm<lb_t, gt_t>::ff_farm(inputCh, inBufferEntries, outBufferEntries, workerCleanup, maxNumWorkers, fixedSize){
    construct();
}

template <typename lb_t, typename gt_t>
adp_ff_farm<lb_t, gt_t>::adp_ff_farm(AdaptivityParameters adaptivityParameters, std::vector<ff_node*>& w,
             ff_node* const emitter, ff_node* const collector, bool inputCh):
    ff_farm<lb_t, gt_t>::ff_farm(w, emitter, collector, inputCh){
    construct(adaptivityParameters);
}

template <typename lb_t, typename gt_t>
adp_ff_farm<lb_t, gt_t>::adp_ff_farm(AdaptivityParameters adaptivityParameters, bool inputCh,
                      int inBufferEntries,
                      int outBufferEntries,
                      bool workerCleanup,
                      int maxNumWorkers,
                      bool fixedSize):
    ff_farm<lb_t, gt_t>::ff_farm(inputCh, inBufferEntries, outBufferEntries, workerCleanup, maxNumWorkers, fixedSize){
    construct(adaptivityParameters);
}

template <typename lb_t, typename gt_t>
adp_ff_farm<lb_t, gt_t>::~adp_ff_farm(){
    ;
}

template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::setAdaptivityParameters(AdaptivityParameters adaptivityParameters){
    _adaptivityParameters = adaptivityParameters;
    uint validationRes = _adaptivityParameters.validate();
    if(validationRes != VALIDATION_OK){
        throw std::runtime_error("AdaptiveFarm: invalid AdaptivityParameters: " + utils::intToString(validationRes));
    }
}

template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::firstRunBefore(){
    if(_firstRun){
        svector<ff_node*> workers = ff_farm<lb_t, gt_t>::getWorkers();
        for(size_t i = 0; i < workers.size(); i++){
            _adaptiveWorkers.push_back(static_cast<adp_ff_node*>(workers[i]));
            _adaptiveWorkers.at(i)->initMammutModules(_adaptivityParameters.mammut);
        }

        _adaptiveEmitter = static_cast<adp_ff_node*>(ff_farm<lb_t, gt_t>::getEmitter());
        if(_adaptiveEmitter){
            _adaptiveEmitter->initMammutModules(_adaptivityParameters.mammut);
        }

        _adaptiveCollector = static_cast<adp_ff_node*>(ff_farm<lb_t, gt_t>::getCollector());
        if(_adaptiveCollector){
            _adaptiveCollector->initMammutModules(_adaptivityParameters.mammut);
        }
    }
}

template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::firstRunAfter(){
    if(_firstRun){
        _firstRun = false;
        _adaptivityManager = new AdaptivityManagerFarm(this, _adaptivityParameters);
        _adaptivityManager->start();
    }
}

template <typename lb_t, typename gt_t>
int adp_ff_farm<lb_t, gt_t>::run(bool skip_init){
    firstRunBefore();
    int r = ff_farm<lb_t, gt_t>::run(skip_init);
    if(r){
        return r;
    }
    firstRunAfter();
    return r;
}

template <typename lb_t, typename gt_t>
int adp_ff_farm<lb_t, gt_t>::wait(){
    int r = 0;
    //int r = ff_farm<lb_t, gt_t>::wait();
    if(_adaptivityManager){
        //_adaptivityManager->stop();
        _adaptivityManager->join();
        delete _adaptivityManager;
    }
    return r;
}


template <typename lb_t, typename gt_t>
void adp_ff_farm<lb_t, gt_t>::waitInternal(){
    ff_farm<lb_t, gt_t>::wait();
}

}

#endif /* ADAPTIVE_FASTFLOW_FARM_HPP_ */
