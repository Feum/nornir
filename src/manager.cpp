/*
 * manager.cpp
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

#include "./manager.hpp"
#include "parameters.hpp"
#include "predictors.hpp"
#include "./node.hpp"
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

class Parameters;
class ManagerFarm;

using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

bool ManagerFarm::reconfigureWorkers() const{
    return true;
}

bool ManagerFarm::reconfigureFrequency() const{
    return _p.knobFrequencies == KNOB_FREQUENCY_YES;
}

uint ManagerFarm::getConfigurationDimension() const{
    uint numDimensions = 0;
    if(reconfigureWorkers()){++numDimensions;}
    if(reconfigureFrequency()){++numDimensions;}
    return numDimensions;
}

vector<PhysicalCore*> ManagerFarm::getSepDomainsPhyCores(
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

void ManagerFarm::setDomainToHighestFrequency(const Domain* domain){
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

vector<VirtualCore*> ManagerFarm::getAvailableVirtualCores(){
    vector<VirtualCore*> r;
    if(_p.knobMapping == KNOB_MAPPING_AUTO){
        _p.knobMapping = KNOB_MAPPING_LINEAR;
    }

    switch(_p.knobMapping){
        case KNOB_MAPPING_LINEAR:{
           /*
            * Generates a vector of virtual cores to be used for linear
            * mapping.node. It contains first one virtual core per physical
            * core (virtual cores on the same CPU are consecutive).
            * Then, the other groups of virtual cores follow.
            */
            vector<Cpu*> cpus = _topology->getCpus();

            size_t virtualUsed = 0;
            size_t virtualPerPhysical;
            if(_p.knobHyperthreading != KNOB_HT_NO){
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
        case KNOB_MAPPING_CACHE_EFFICIENT:{
            throw runtime_error("Not yet supported.");
        }
        default:
            break;
    }
    return r;
}

uint ManagerFarm::numScalableServiceNodes(){
    return (_emitter &&
            _p.knobMappingEmitter != KNOB_SNODE_MAPPING_PERFORMANCE) +
           (_collector &&
            _p.knobMappingEmitter != KNOB_SNODE_MAPPING_PERFORMANCE);
}

void ManagerFarm::manageServiceNodesPerformance(){
    if((_p.knobMappingEmitter != KNOB_SNODE_MAPPING_ALONE) && !_emitter){
        _p.knobMappingEmitter = KNOB_SNODE_MAPPING_ALONE;
    }

    if((_p.knobMappingCollector != KNOB_SNODE_MAPPING_ALONE) && !_collector){
        _p.knobMappingCollector = KNOB_SNODE_MAPPING_ALONE;
    }

    if((_p.knobMappingEmitter == KNOB_SNODE_MAPPING_PERFORMANCE &&
            !_emitterSensitivitySatisfied) ||
       (_p.knobMappingCollector == KNOB_SNODE_MAPPING_PERFORMANCE &&
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

            if(_p.knobMappingEmitter == KNOB_SNODE_MAPPING_PERFORMANCE){
                VirtualCore* vc = perfPhyCores.at(index)->getVirtualCore();
                setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                _emitterVirtualCore = vc;
                _emitterSensitivitySatisfied = true;
                index = (index + 1) % perfPhyCores.size();
            }

            if(_p.knobMappingCollector == KNOB_SNODE_MAPPING_PERFORMANCE){
                VirtualCore* vc = perfPhyCores.at(index)->getVirtualCore();
                setDomainToHighestFrequency(_cpufreq->getDomain(vc));
                _collectorVirtualCore = vc;
                _collectorSensitivitySatisfied = true;
                index = (index + 1) % perfPhyCores.size();
            }
        }
    }
}

void ManagerFarm::getMappingIndexes(size_t& emitterIndex,
                       size_t& firstWorkerIndex,
                       size_t& collectorIndex){
    size_t nextIndex = 0;
    if(_emitter && !_emitterVirtualCore){
        emitterIndex = nextIndex;
        if(_p.knobMappingEmitter != KNOB_SNODE_MAPPING_COLLAPSED){
            nextIndex = (nextIndex + 1) % _availableVirtualCores.size();
        }
    }

    firstWorkerIndex = nextIndex;
    nextIndex = (nextIndex + _activeWorkers.size()) %
                _availableVirtualCores.size();

    if(_collector && !_collectorVirtualCore){
        if(_p.knobMappingCollector == KNOB_SNODE_MAPPING_COLLAPSED){
            nextIndex = (nextIndex - 1) % _availableVirtualCores.size();
        }
        collectorIndex = nextIndex;
    }
}

void ManagerFarm::mapNodesToVirtualCores(){
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

void ManagerFarm::applyUnusedVCStrategyOff(const vector<VirtualCore*>& unusedVc){
    for(size_t i = 0; i < unusedVc.size(); i++){
        VirtualCore* vc = unusedVc.at(i);
        if(vc->isHotPluggable() && vc->isHotPlugged()){
            vc->hotUnplug();
        }
    }
}

void ManagerFarm::applyUnusedVCStrategyLowestFreq(const vector<VirtualCore*>& vc){
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

void ManagerFarm::applyUnusedVCStrategy(StrategyUnusedVirtualCores strategyUnused,
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

void ManagerFarm::applyUnusedVCStrategy(){
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

void ManagerFarm::updateScalableDomains(){
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

void ManagerFarm::updatePstate(Frequency frequency){
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
        }
    }
}

void ManagerFarm::mapAndSetFrequencies(){
    if(_p.knobMapping == KNOB_MAPPING_NO){
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

    if(reconfigureFrequency()){
        // We suppose that all the domains have the same
        // available frequencies.
        _availableFrequencies = _cpufreq->getDomains().at(0)->
                                getAvailableFrequencies();

        // Remove turbo boost frequency.
        if(!_p.turboBoost && !_availableFrequencies.empty()){
            if(intToString(_availableFrequencies.back()).at(3) == '1'){
                _availableFrequencies.pop_back();
            }
        }
    }

    if(_availableFrequencies.empty()){
        // Insert dummy constant frequency
        _availableFrequencies.push_back(1.0);
    }

    // Sets the current frequency to the highest possible.
    _currentConfiguration.frequency = _availableFrequencies.back();

    if(_p.knobFrequencies != KNOB_FREQUENCY_NO){
        updatePstate(_currentConfiguration.frequency);
    }
}

double ManagerFarm::getMaxPredictionErrorPrimary() const{
    double r = _p.maxPrimaryPredictionError;
    if(_p.strategyPredictionErrorPrimary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getPrimaryValue(_samples->coefficientVariation()));
    }
    return r;
}

double ManagerFarm::getMaxPredictionErrorSecondary() const{
    double r = _p.maxSecondaryPredictionError;
    if(_p.strategyPredictionErrorSecondary ==
       STRATEGY_PREDICTION_ERROR_COEFFVAR){
        r = max(r, getSecondaryValue(_samples->coefficientVariation()));
    }
    return r;
}

double ManagerFarm::getPrimaryValue(const MonitoredSample& sample) const{
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

double ManagerFarm::getSecondaryValue(const MonitoredSample& sample) const{
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

double ManagerFarm::getPrimaryValue() const{
    return getPrimaryValue(_samples->average());
}

double ManagerFarm::getSecondaryValue() const{
    return getSecondaryValue(_samples->average());
}

bool ManagerFarm::isContractViolated() const{
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

bool ManagerFarm::isFeasiblePrimaryValue(double value, double tolerance) const{
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

double ManagerFarm::getVoltage(const FarmConfiguration& configuration) const{
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

bool ManagerFarm::isBestSuboptimalValue(double x, double y) const{
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

bool ManagerFarm::isBestSecondaryValue(double x, double y) const{
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

FarmConfiguration ManagerFarm::getNewConfiguration(){
    FarmConfiguration bestConfiguration;
    FarmConfiguration bestSuboptimalConfiguration = _currentConfiguration;

    double primaryPrediction = 0;
    double secondaryPrediction = 0;

    double bestPrimaryPrediction = 0;
    double bestSecondaryPrediction = 0;
    double bestSuboptimalValue = getPrimaryValue();

    bool feasibleSolutionFound = false;

    switch(_p.contractType){
        case CONTRACT_PERF_UTILIZATION:
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:{
            // We have to minimize the power/energy.
            bestSecondaryPrediction = numeric_limits<double>::max();
        }break;
        case CONTRACT_POWER_BUDGET:{
            // We have to maximize the bandwidth.
            bestSecondaryPrediction = numeric_limits<double>::min();
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
                                        bestSecondaryPrediction)){
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
            }
        }
    }

    if(feasibleSolutionFound){
        _primaryPrediction = bestPrimaryPrediction;
        _secondaryPrediction = bestSecondaryPrediction;
        return bestConfiguration;
    }else{
        _primaryPrediction = bestSuboptimalValue;
        // TODO: This check now works because both service time and power are always  > 0
        // In the future we must find another way to indicat that secondary prediction
        // has not been done.
        _secondaryPrediction = -1;
        return bestSuboptimalConfiguration;
    }
}

void ManagerFarm::updateUsedCpus(){
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

void ManagerFarm::changeActiveNodes(FarmConfiguration& configuration){
    uint workersNumDiff = abs((int)_currentConfiguration.numWorkers -
                              (int)configuration.numWorkers);
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

void ManagerFarm::prepareToFreeze(){
    if(_emitter){
        _emitter->prepareToFreeze();
    }
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->prepareToFreeze();
    }
    if(_collector){
        _collector->prepareToFreeze();
    }
}

void ManagerFarm::freeze(){
    if(_emitter){
        _emitter->freezeAll((void*)(_collector?FF_EOSW:FF_GO_OUT));
    }
    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->wait_freezing();
    }
    if(_collector){
        _collector->wait_freezing();
    }
}

void ManagerFarm::notifyNewConfiguration(uint numWorkers){
    if(_emitter){
        _emitter->notifyWorkersChange(_currentConfiguration.numWorkers,
                                      numWorkers);
    }
    for(size_t i = 0; i < numWorkers; i++){
        AdaptiveNode* w = _activeWorkers.at(i);
        w->notifyWorkersChange(_currentConfiguration.numWorkers,
                               numWorkers);
    }
    if(_collector){
        _collector->notifyWorkersChange(_currentConfiguration.numWorkers,
                                        numWorkers);
    }
}

void ManagerFarm::prepareToRun(uint numWorkers){
    if(_emitter){
        _emitter->prepareToRun();
    }
    for(size_t i = 0; i < numWorkers; i++){
        _activeWorkers.at(i)->prepareToRun();
    }
    if(_collector){
        _collector->prepareToRun();
    }
}

void ManagerFarm::run(uint numWorkers){
    if(_emitter){
        _emitter->thawAll(numWorkers);
    }
    if(_collector){
        _farm->getgt()->thaw(true, numWorkers);
    }
}

bool ManagerFarm::terminated(){
    if(_emitter &&
       _emitter->isTerminated()){
        return true;
    }

    for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
        if(_activeWorkers.at(i)->isTerminated()){
            return true;
        }
    }

    if(_collector &&
       _collector->isTerminated()){
        return true;
    }

    return false;
}

void ManagerFarm::changeConfiguration(FarmConfiguration configuration){
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

        prepareToFreeze();
        freeze();

        changeActiveNodes(configuration);
        updateUsedCpus();
        notifyNewConfiguration(configuration.numWorkers);

        prepareToRun(configuration.numWorkers);
        run(configuration.numWorkers);

        if(_p.fastReconfiguration){
            _cpufreq->rollback(rollbackPoints);
        }

        // Must be done after rollback.
        applyUnusedVCStrategy();
    }
    /****************** Workers change terminated ******************/


    /****************** P-state change started ******************/
    //TODO: Maybe sensitivity could not be satisfied with the maximum
    // number of workers but could be satisfied now.
    if(_p.knobFrequencies != KNOB_FREQUENCY_NO){
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

void ManagerFarm::observe(){
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

void ManagerFarm::askForWorkersSamples(){
    for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
        _activeWorkers.at(i)->askForSample();
    }
}

void ManagerFarm::getWorkersSamples(WorkerSample& sample){
    AdaptiveNode* w;
    sample = WorkerSample();
    for(size_t i = 0; i < _currentConfiguration.numWorkers; i++){
        WorkerSample tmp;
        w = _activeWorkers.at(i);
        w->getSampleResponse(tmp, _p.strategyPolling,
                             _samples->average().latency);
        sample += tmp;
    }
    sample.loadPercentage /= _currentConfiguration.numWorkers;
    sample.latency /= _currentConfiguration.numWorkers;
}

void ManagerFarm::storeNewSample(){
    MonitoredSample sample;
    WorkerSample ws;
    JoulesCpu joules;

    askForWorkersSamples();
    getWorkersSamples(ws);

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
}

bool ManagerFarm::persist() const{
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

void ManagerFarm::initCalibrator(){
    if(_p.strategyCalibration == STRATEGY_CALIBRATION_RANDOM){
        ;
    }else{
        _calibrator = new CalibratorLowDiscrepancy(*this);
    }
}

void ManagerFarm::initPredictors(){
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

ManagerFarm::ManagerFarm(ff_farm<>* farm, Parameters parameters):
        _farm(farm),
        _p(parameters),
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

    ParametersValidation apv = _p.validate();
    if(apv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + apv);
    }
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

ManagerFarm::~ManagerFarm(){
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

void ManagerFarm::run(){
    _emitter = static_cast<AdaptiveNode*>(_farm->getEmitter());
    _collector = static_cast<AdaptiveNode*>(_farm->getCollector());
    svector<ff_node*> w = _farm->getWorkers();
    for(size_t i = 0; i < w.size(); i++){
        _activeWorkers.push_back(static_cast<AdaptiveNode*>(w[i]));
    }
    _maxNumWorkers = _activeWorkers.size();
    _currentConfiguration = FarmConfiguration(_maxNumWorkers, 0);

    _farm->run_then_freeze(_maxNumWorkers);

    for(size_t i = 0; i < _activeWorkers.size(); i++){
        _activeWorkers.at(i)->init(_p.mammut, _p.archData.ticksPerNs);
    }
    if(_emitter){
        _emitter->init(_p.mammut, _p.archData.ticksPerNs);
    }else{
        throw runtime_error("Emitter is needed to use the manager.");
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
            changeConfiguration(_calibrator->getNextConfiguration());
        }

        double startSample = getMillisecondsTime();
        while(!terminated()){
            double overheadMs = getMillisecondsTime() - startSample;
            microsecsSleep = ((double)_p.samplingInterval - overheadMs)*
                              (double)MAMMUT_MICROSECS_IN_MILLISEC;
            if(microsecsSleep < 0){
                microsecsSleep = 0;
            }
            usleep(microsecsSleep);
            startSample = getMillisecondsTime();

            storeNewSample();

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

double Observer::calibrationDurationToPerc(const CalibrationStats& cs,
                                 uint durationMs){
    return ((double)cs.duration /
            (double)durationMs) * 100.0;
}

Observer::Observer(string statsFile, string calibrationFile, string summaryFile):
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

Observer::~Observer(){
    _statsFile.close();
    _calibrationFile.close();
    _summaryFile.close();
}

void Observer::observe(unsigned int timeStamp,
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

void Observer::calibrationStats(const vector<CalibrationStats>&
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

void Observer::summaryStats(const vector<CalibrationStats>&
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


}
