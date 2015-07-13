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
 *  1. Emitter, Workers and Collector of the farm must extend adpff::adp_ff_node
 *     instead of ff::ff_node
 *  2. Replace the following calls (if present) in the farm nodes:
 *          svc           -> adp_svc
 *          svc_init      -> adp_svc_init
 *          svc_end       -> adp_svc_end
 *  3. If the application wants to be aware of the changes in the number of workers, the nodes
 *     can implement the notifyWorkersChange virtual method.
 *  4. Substitute ff::ff_farm with adpff::adp_ff_farm. The maximum number of workers
 *     that can be activated correspond to the number of workers specified during farm creation.
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

class AdaptivityParameters;
class adp_ff_node;

template<typename lb_t, typename gt_t>
class AdaptivityManagerFarm;

using namespace ff;
using namespace mammut;

/*!
 * This class can be used to obtain statistics about reconfigurations
 * performed by the manager.
 * It can be extended by a user defined class to customize action to take
 * for each observed statistic.
 */
class adp_ff_farm_observer{
    template<typename lb_t, typename gt_t>
    friend class AdaptivityManagerFarm;
protected:
    size_t _numberOfWorkers;
    cpufreq::Frequency _currentFrequency;
    topology::VirtualCore* _emitterVirtualCore;
    std::vector<topology::VirtualCore*> _workersVirtualCore;
    topology::VirtualCore* _collectorVirtualCore;
    double _currentTasks;
    double _currentUtilization;
    double _currentCoresUtilization;
    energy::JoulesCpu _usedJoules;
    energy::JoulesCpu _unusedJoules;
    unsigned int _startMonitoringMs;
public:
    adp_ff_farm_observer():_numberOfWorkers(0), _currentFrequency(0), _emitterVirtualCore(NULL),
                         _collectorVirtualCore(NULL), _currentTasks(0), _currentUtilization(0),
                         _currentCoresUtilization(0){;}
    virtual ~adp_ff_farm_observer(){;}
    virtual void observe(){;}
};

/*!
 * \class adp_ff_farm
 * \brief This class wraps a farm to let it reconfigurable.
 *
 * This class wraps a farm to let it reconfigurable.
 */
template<typename lb_t=ff_loadbalancer, typename gt_t=ff_gatherer>
class adp_ff_farm: public ff_farm<lb_t, gt_t>{
    friend class AdaptivityManagerFarm<lb_t, gt_t>;
public:
    void setAdaptivityParameters(AdaptivityParameters adaptivityParameters){
		_adaptivityParameters = adaptivityParameters;
		uint validationRes = _adaptivityParameters.validate();
		if(validationRes != VALIDATION_OK){
			throw std::runtime_error("AdaptiveFarm: invalid AdaptivityParameters: " + utils::intToString(validationRes));
		}
	}
private:
    std::vector<adp_ff_node*> _adaptiveWorkers;
    adp_ff_node* _adaptiveEmitter;
    adp_ff_node* _adaptiveCollector;
    bool _firstRun;
    AdaptivityParameters _adaptivityParameters;
    AdaptivityManagerFarm<lb_t, gt_t>* _adaptivityManager;

    void construct(){
        _adaptiveEmitter = NULL;
        _adaptiveCollector = NULL;
        _firstRun = true;
        _adaptivityManager = NULL;
    }

    void construct(AdaptivityParameters adaptivityParameters){
    	construct();
        setAdaptivityParameters(adaptivityParameters);
    }

    std::vector<adp_ff_node*> getAdaptiveWorkers() const{
        return _adaptiveWorkers;
    }

    adp_ff_node* getAdaptiveEmitter() const{
        return _adaptiveEmitter;
    }

    adp_ff_node* getAdaptiveCollector() const{
        return _adaptiveCollector;
    }
public:

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     */
    adp_ff_farm(std::vector<ff_node*>& w, ff_node* const emitter = NULL,
    		    ff_node* const collector = NULL, bool inputCh = false):
		ff_farm<lb_t, gt_t>::ff_farm(w, emitter, collector, inputCh){
		construct();
    }

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     */
    adp_ff_farm(bool inputCh = false,
                          int inBufferEntries = ff_farm<lb_t, gt_t>::DEF_IN_BUFF_ENTRIES,
                          int outBufferEntries = ff_farm<lb_t, gt_t>::DEF_OUT_BUFF_ENTRIES,
                          bool workerCleanup = false,
                          int maxNumWorkers = ff_farm<lb_t, gt_t>::DEF_MAX_NUM_WORKERS,
                          bool fixedSize = false):
		ff_farm<lb_t, gt_t>::ff_farm(inputCh, inBufferEntries, outBufferEntries, workerCleanup, maxNumWorkers, fixedSize){
		construct();
	}

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     * @param adaptivityParameters Parameters that will be used by the farm to take reconfiguration decisions.
     */
    adp_ff_farm(AdaptivityParameters adaptivityParameters, std::vector<ff_node*>& w,
                 ff_node* const emitter = NULL, ff_node* const collector = NULL, bool inputCh = false):
		ff_farm<lb_t, gt_t>::ff_farm(w, emitter, collector, inputCh){
		construct(adaptivityParameters);
	}

    /**
     * Builds the adaptive farm.
     * For parameters documentation, see fastflow's farm documentation.
     * @param adaptivityParameters Parameters that will be used by the farm to take reconfiguration decisions.
     */
    explicit adp_ff_farm(AdaptivityParameters adaptivityParameters, bool inputCh = false,
                          int inBufferEntries = ff_farm<lb_t, gt_t>::DEF_IN_BUFF_ENTRIES,
                          int outBufferEntries = ff_farm<lb_t, gt_t>::DEF_OUT_BUFF_ENTRIES,
                          bool workerCleanup = false,
                          int maxNumWorkers = ff_farm<lb_t, gt_t>::DEF_MAX_NUM_WORKERS,
                          bool fixedSize = false):
		ff_farm<lb_t, gt_t>::ff_farm(inputCh, inBufferEntries, outBufferEntries, workerCleanup, maxNumWorkers, fixedSize){
		construct(adaptivityParameters);
	}

    /**
     * Destroyes this adaptive farm.
     */
    ~adp_ff_farm(){
        ;
    }

    void firstRunBefore(){
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


    void firstRunAfter(){
        if(_firstRun){
            _firstRun = false;
            _adaptivityManager = new AdaptivityManagerFarm<lb_t, gt_t>(this, _adaptivityParameters);
            _adaptivityManager->start();
        }
    }

    /**
     * Runs this farm.
     */
    int run(bool skip_init=false){
    	firstRunBefore();
        int r = ff_farm<lb_t, gt_t>::run(skip_init);
        if(r){
            return r;
        }
        firstRunAfter();
        return r;
    }

    /**
     * Waits this farm for completion.
     */
    int wait(){
        int r = ff_farm<lb_t, gt_t>::wait();
        if(_adaptivityManager){
            _adaptivityManager->stop();
            _adaptivityManager->join();
            delete _adaptivityManager;
        }
        return r;
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
    cpufreq::Frequency frequency;
    FarmConfiguration():numWorkers(0), frequency(0){;}
    FarmConfiguration(uint numWorkers, cpufreq::Frequency frequency = 0):numWorkers(numWorkers), frequency(frequency){;}
}FarmConfiguration;

typedef struct MonitoredSample{
    std::vector<NodeSample> nodes; ///< The samples taken from the active workers (one per worker).
    energy::JoulesCpu usedCpusEnergy; ///< Energy consumed by active CPUs.
    energy::JoulesCpu unusedCpusEnergy; ///< Energy consumed by inactive CPUs.


    NodeSample getAverageNodesSamples() const{
        NodeSample ns;
        for(size_t i = 0; i < nodes.size(); i++){
            ns += nodes.at(i);
        }
        ns /= nodes.size();
        return ns;
    }

    MonitoredSample& operator+=(const MonitoredSample& rhs){
        for(size_t i = 0; i < nodes.size(); i++){
            nodes.at(i) = nodes.at(i) + rhs.nodes.at(i);
        }
        usedCpusEnergy += rhs.usedCpusEnergy;
        unusedCpusEnergy += rhs.unusedCpusEnergy;
        // actual addition of rhs to *this
        return *this;
    }

    MonitoredSample& operator/=(double c){
        for(size_t i = 0; i < nodes.size(); i++){
            nodes.at(i) = nodes.at(i) / c;
        }
        usedCpusEnergy /= c;
        unusedCpusEnergy /= c;
        return *this;

    }
}MonitoredSample;

inline MonitoredSample operator+(MonitoredSample lhs, const MonitoredSample& rhs){
    lhs += rhs;
    return lhs;
}

inline MonitoredSample operator/(MonitoredSample lhs, double c){
    lhs /= c;
    return lhs;
}

/*!
 * \internal
 * \class AdaptivityManagerFarm
 * \brief This class manages the adaptivity in farm based computations.
 *
 * This class manages the adaptivity in farm based computations.
 */
template<typename lb_t=ff_loadbalancer, typename gt_t=ff_gatherer>
class AdaptivityManagerFarm: public utils::Thread{
    friend class PredictorSimple;
    friend class PredictorLinearRegression;
private:
    utils::Monitor _monitor; ///< Used to let the manager stop safe.
    adp_ff_farm<lb_t, gt_t>* _farm; ///< The managed farm.
    AdaptivityParameters _p; ///< The parameters used to take management decisions.
    cpufreq::CpuFreq* _cpufreq; ///< The cpufreq module.
    energy::Energy* _energy; ///< The energy module.
    task::TasksManager* _task; ///< The task module.
    topology::Topology* _topology; ///< The topology module.
    uint _numPhysicalCores; ///< Number of physical cores.
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
    Window<MonitoredSample> _monitoredSamples; ///< Monitored samples;
    double _averageTasks; ///< The last value registered for average tasks per sampling interval processed.
    double _averageUtilization; ///< The last value registered for average utilization per sampling interval.
    double _averageCoresUtilization; ///< The last value registered for average cores utilization (percentage of the average utilization).
    energy::JoulesCpu _usedJoules; ///< Average Joules consumed by the used virtual cores per sampling interval.
    energy::JoulesCpu _unusedJoules; ///< Average Joules consumed by the unused virtual cores per sampling interval.
    uint64_t _remainingTasks; ///< When contract is CONTRACT_COMPLETION_TIME, represent the number of tasks that
                             ///< still needs to be processed by the application.
    time_t _deadline; ///< When contract is CONTRACT_COMPLETION_TIME, represent the deadline of the application.

    Predictor* _primaryPredictor; ///< The predictor of the primary value.
    Predictor* _secondaryPredictor; ///< The predictor of the secondary value.

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
                size_t virtualPerPhysical = _topology->getVirtualCores().size() /
                                            _topology->getPhysicalCores().size();
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
     * Generates mapping indexes. They are indexes to be used on _availableVirtualCores vector
     * to get the corresponding virtual core where a specific node must be mapped.
     * @param emitterIndex The index of the emitter.
     * @param firstWorkerIndex The index of the first worker (the others follow).
     * @param collectorIndex The index of the collector (if present).
     */
    void getMappingIndexes(size_t& emitterIndex, size_t& firstWorkerIndex, size_t& collectorIndex){
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
         * TODO: Explain better
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
        for(size_t i = 0; i < _scalableDomains.size(); i++){
            if(!_scalableDomains.at(i)->setGovernor(_p.frequencyGovernor)){
                throw std::runtime_error("AdaptivityManagerFarm: Impossible to set the specified governor.");
            }
            if(_p.frequencyGovernor != cpufreq::GOVERNOR_USERSPACE){
                if(!_scalableDomains.at(i)->setGovernorBounds(_p.frequencyLowerBound,
                                                              _p.frequencyUpperBound)){
                    throw std::runtime_error("AdaptivityManagerFarm: Impossible to set the specified governor's bounds.");
                }
            }else if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                if(!_scalableDomains.at(i)->setFrequencyUserspace(frequency)){
                    throw std::runtime_error("AdaptivityManagerFarm: Impossible to set the specified frequency.");
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

        _availableFrequencies.push_back(1.0); // Insert dummy constant frequency
        if(_p.strategyFrequencies != STRATEGY_FREQUENCY_NO){
            if(_p.strategyFrequencies != STRATEGY_FREQUENCY_OS){
                /** We suppose that all the domains have the same available frequencies. **/
                _availableFrequencies = _cpufreq->getDomains().at(0)->getAvailableFrequencies();
                /** Sets the current frequency to the highest possible. **/
                _currentConfiguration.frequency = _availableFrequencies.back();
            }
            updatePstate(_currentConfiguration.frequency);
        }
    }

    /**
     * Updates the monitored values.
     */
    void updateMonitoredValues(){
        _averageTasks = 0;
        _averageUtilization = 0;
        _averageCoresUtilization = 0;
        _usedJoules.zero();
        _unusedJoules.zero();

        MonitoredSample average = _monitoredSamples.average();
        NodeSample nodesAvg = average.getAverageNodesSamples();
        _averageTasks = nodesAvg.tasksCount * _currentConfiguration.numWorkers;
        _averageCoresUtilization = nodesAvg.corePercentage;
        _averageUtilization = nodesAvg.loadPercentage;
        _usedJoules = average.usedCpusEnergy;
        _unusedJoules = average.unusedCpusEnergy;
    }

    /**
     * Returns the monitored value according to the required contract.
     * @return The monitored value according to the required contract.
     */
    double getMonitoredValue() const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return _averageUtilization;
            }break;
            case CONTRACT_PERF_BANDWIDTH:
            case CONTRACT_PERF_COMPLETION_TIME:{
                return _averageTasks / (double) _p.samplingInterval;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return (_usedJoules.cores + _unusedJoules.cores) / (double) _p.samplingInterval;
            }break;
            default:{
                return 0;
            }break;
        }
    }

    /**
     * Checks if a specific monitored value violates the contract requested
     * by the user.
     * @param monitoredValue The monitored value.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated(double monitoredValue) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                return monitoredValue < _p.underloadThresholdFarm ||
                       monitoredValue > _p.overloadThresholdFarm;
            }break;
            case CONTRACT_PERF_BANDWIDTH:{
                double offset = (_p.requiredBandwidth * _p.maxBandwidthVariation) / 100.0;
                return monitoredValue < _p.requiredBandwidth ||
                       monitoredValue > _p.requiredBandwidth + offset;
            }break;
            case CONTRACT_PERF_COMPLETION_TIME:{
                return monitoredValue < _p.requiredBandwidth;
            }break;
            case CONTRACT_POWER_BUDGET:{
                return monitoredValue < _p.powerBudget;
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
        cpufreq::VoltageTableKey key(configuration.numWorkers, configuration.frequency);
        cpufreq::VoltageTableIterator it = _voltageTable.find(key);
        if(it != _voltageTable.end()){
            return it->second;
        }else{
            throw std::runtime_error("Frequency and/or number of virtual cores not found in voltage table.");
        }
    }

    /**
     * Checks if x is a best suboptimal monitored value than y.
     * @param x The first monitored value.
     * @param y The second monitored value.
     * @return True if x is a best suboptimal monitored value than y, false otherwise.
     */
    bool isBestSuboptimalValue(double x, double y) const{
        switch(_p.contractType){
            case CONTRACT_PERF_UTILIZATION:{
                // Concerning utilization factors, if both are suboptimal, we prefer the closest to the lower bound.
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
                // Concerning bandwidths, if both are suboptimal, we prefer the higher one.
                return x > y;
            }break;
            case CONTRACT_POWER_BUDGET:{
                // Concerning power budgets, if both are suboptimal, we prefer the lowest one.
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
    FarmConfiguration getNewConfiguration() const{
        FarmConfiguration r;
        double predictedMonitoredValue = 0;
        double bestSuboptimalValue = getMonitoredValue();
        FarmConfiguration bestSuboptimalConfiguration = _currentConfiguration;
        bool feasibleSolutionFound = false;

        double bestSecondaryValue;
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
        double predictedSecondaryValue = 0;
        unsigned int remainingTime = 0;
        for(size_t i = 1; i <= _maxNumWorkers; i++){
            for(size_t j = 0; j < _availableFrequencies.size(); j++){
                FarmConfiguration examinedConfiguration(i, _availableFrequencies.at(j));
                predictedMonitoredValue = _primaryPredictor->predict(examinedConfiguration);
                switch(_p.contractType){
                    case CONTRACT_PERF_COMPLETION_TIME:{
                        remainingTime = (double) _remainingTasks / predictedMonitoredValue;
                    }break;
                    case CONTRACT_PERF_UTILIZATION:{
                        predictedMonitoredValue = ((_averageTasks / _p.samplingInterval) / predictedMonitoredValue) * _averageUtilization;
                    }break;
                    default:{
                        ;
                    }
                }

                if(!isContractViolated(predictedMonitoredValue)){
                    predictedSecondaryValue = _secondaryPredictor->predict(examinedConfiguration);
                    if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                        predictedSecondaryValue *= remainingTime;
                    }
                    if(isBestSecondaryValue(predictedSecondaryValue, bestSecondaryValue)){
                        bestSecondaryValue = predictedSecondaryValue;
                        r = examinedConfiguration;
                        feasibleSolutionFound = true;
                    }
                }else if(!feasibleSolutionFound && isBestSuboptimalValue(predictedMonitoredValue, bestSuboptimalValue)){
                    bestSuboptimalValue = predictedMonitoredValue;
                    bestSuboptimalConfiguration = examinedConfiguration;
                }
            }
        }

        if(feasibleSolutionFound){
            return r;
        }else{
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
     */
    void changeConfiguration(FarmConfiguration configuration){
        if(!configuration.numWorkers){
            throw std::runtime_error("AdaptivityManagerFarm: fatal error, trying to activate zero "
                                     "workers.");
        }else if(configuration.numWorkers > _maxNumWorkers){
            throw std::runtime_error("AdaptivityManagerFarm: fatal error, trying to activate more "
                                     "workers than the maximum allowed.");
        }

        /****************** Refine the model ******************/
        _primaryPredictor->refine();
        _secondaryPredictor->refine();

        /****************** Workers change started ******************/
        if(_currentConfiguration.numWorkers != configuration.numWorkers){
            std::vector<cpufreq::RollbackPoint> rollbackPoints;
            if(_p.fastReconfiguration){
                for(size_t i = 0; i < _scalableDomains.size(); i++){
                    rollbackPoints.push_back(_scalableDomains.at(i)->getRollbackPoint());
                    setDomainToHighestFrequency(_scalableDomains.at(i));
                }
            }


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

            /** Stops farm. **/
            _emitter->produceNull();
            _farm->wait_freezing();
            /** Notify the nodes that a reconfiguration is happening. **/
            _emitter->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            for(uint i = 0; i < configuration.numWorkers; i++){
                _activeWorkers.at(i)->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            }
            if(_collector){
                _collector->notifyWorkersChange(_currentConfiguration.numWorkers, configuration.numWorkers);
            }
            /** Start the farm again. **/
            _farm->run_then_freeze(configuration.numWorkers);
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
    }

    /** Send data to observer. **/
    void observe(){
        if(_p.observer){
            _p.observer->_numberOfWorkers = _currentConfiguration.numWorkers;
            _p.observer->_currentFrequency = _currentConfiguration.frequency;
            _p.observer->_emitterVirtualCore = _emitterVirtualCore;
            _p.observer->_workersVirtualCore = _activeWorkersVirtualCores;
            _p.observer->_collectorVirtualCore = _collectorVirtualCore;
            _p.observer->_currentTasks = _averageTasks;
            _p.observer->_currentUtilization = _averageUtilization;
            _p.observer->_currentCoresUtilization = _averageCoresUtilization;
            _p.observer->_usedJoules = _usedJoules;
            _p.observer->_unusedJoules = _unusedJoules;
            _p.observer->observe();
        }
    }

    /**
     * Asks the workers for samples.
     * @param nextSampleIndex The index where to store the new sample.
     * @return false if the farm is not running anymore, true otherwise.
     **/
    bool storeNewSamples(){
        for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
            _activeWorkers.at(i)->askForSample();
        }

        MonitoredSample sample;
        for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
            NodeSample ns;
            bool workerRunning = _activeWorkers.at(i)->getSampleResponse(ns);
            if(!workerRunning){
                return false;
            }

            if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
                if(_remainingTasks > ns.tasksCount){
                    _remainingTasks -= ns.tasksCount;
                }else{
                    _remainingTasks = 0;
                }
            }

            sample.nodes.push_back(ns);
        }

        sample.usedCpusEnergy.zero();
        for(size_t i = 0; i < _usedCpus.size(); i++){
            energy::CounterCpu* currentCounter = _energy->getCounterCpu(_usedCpus.at(i));
            sample.usedCpusEnergy += currentCounter->getJoules();
        }
        sample.unusedCpusEnergy.zero();
        for(size_t i = 0; i < _unusedCpus.size(); i++){
            energy::CounterCpu* currentCounter = _energy->getCounterCpu(_unusedCpus.at(i));
            sample.unusedCpusEnergy += currentCounter->getJoules();
        }
        _energy->resetCountersCpu();

        _monitoredSamples.add(sample);

        return true;
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
            }break;
            case STRATEGY_PREDICTION_REGRESSION_LINEAR:{
                _primaryPredictor = new PredictorLinearRegression(primary, *this);
                _secondaryPredictor = new PredictorLinearRegression(secondary, *this);
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
    AdaptivityManagerFarm(adp_ff_farm<lb_t, gt_t>* farm, AdaptivityParameters adaptivityParameters):
        _farm(farm),
        _p(adaptivityParameters),
        _cpufreq(_p.mammut.getInstanceCpuFreq()),
        _energy(_p.mammut.getInstanceEnergy()),
        _task(_p.mammut.getInstanceTask()),
        _topology(_p.mammut.getInstanceTopology()),
        _numPhysicalCores(_topology->getPhysicalCores().size()),
        _emitter(_farm->getAdaptiveEmitter()),
        _collector(_farm->getAdaptiveCollector()),
        _activeWorkers(_farm->getAdaptiveWorkers()),
        _maxNumWorkers(_activeWorkers.size()),
        _emitterSensitivitySatisfied(false),
        _collectorSensitivitySatisfied(false),
        _currentConfiguration(_activeWorkers.size(), 0),
        _availableVirtualCores(getAvailableVirtualCores()),
        _emitterVirtualCore(NULL),
        _collectorVirtualCore(NULL),
        _monitoredSamples(_p.numSamples){
        /** If voltage table file is specified, then load the table. **/
        if(_p.voltageTableFile.compare("")){
            cpufreq::loadVoltageTable(_voltageTable, _p.voltageTableFile);
        }
    }

    /**
     * Destroyes this adaptivity manager.
     */
    ~AdaptivityManagerFarm(){
        ;
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

        if(_cpufreq->isBoostingSupported()){
            if(_p.turboBoost){
                _cpufreq->enableBoosting();
            }else{
                _cpufreq->disableBoosting();
            }
        }

        mapAndSetFrequencies();

        _energy->resetCountersCpu();
        _p.observer->_startMonitoringMs = utils::getMillisecondsTime();

        uint64_t samplesToDiscard = _p.samplesToDiscard;
        double lastOverheadMs = 0;
        double startOverheadMs = 0;
        double microsecsSleep = 0;
        if(_p.contractType == CONTRACT_PERF_COMPLETION_TIME){
            _remainingTasks = _p.expectedTasksNumber;
            _deadline = time(NULL) + _p.requiredCompletionTime;
        }

        initPredictors();

        if(_p.contractType == CONTRACT_NONE){
            _monitor.wait();
            storeNewSamples();
            updateMonitoredValues();
            observe();
        }else{
            while(!mustStop()){
                microsecsSleep = (double)_p.samplingInterval*(double)MAMMUT_MICROSECS_IN_SEC -
                                 lastOverheadMs*(double)MAMMUT_MICROSECS_IN_MILLISEC;
                if(microsecsSleep < 0){
                    microsecsSleep = 0;
                }
                usleep(microsecsSleep);

                startOverheadMs = utils::getMillisecondsTime();

                if(!storeNewSamples()){
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

                if(!samplesToDiscard){
                    updateMonitoredValues();
                    observe();
                }else{
                    --samplesToDiscard;
                }

                if((_monitoredSamples.size() >= _p.numSamples) &&
                   isContractViolated(getMonitoredValue())){
                    changeConfiguration(getNewConfiguration());
                    _monitoredSamples.reset();
                    samplesToDiscard = _p.samplesToDiscard;
                }

                lastOverheadMs = utils::getMillisecondsTime() - startOverheadMs;
            }
        }
    controlLoopEnd:
        ;
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

}

#endif /* ADAPTIVE_FASTFLOW_FARM_HPP_ */
