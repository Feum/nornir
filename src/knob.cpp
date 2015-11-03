/*
 * knob.hpp
 *
 * Created on: 02/11/2015
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

#include "knob.hpp"
#include "parameters.hpp"

#include <mammut/mammut.hpp>
#include <mammut/utils.hpp>
#include <mammut/cpufreq/cpufreq.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace adpff{

using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

using namespace ff;

using namespace std;

void Knob::setRelativeValue(double v){
    // Maps from the range [0, 100] to the real range.
    vector<double> values = getAllowedValues();
    if(values.size()){
        uint index = round((double)(values.size() - 1) *
                           (v / 100.0));

        _relativeValue = v;
        setRealValue(values.at(index));
    }
}

void Knob::setRealValue(double v){
    if(getAllowedValues().size() &&
       v != _realValue){
        _realValue = v;
        changeValueReal(_realValue);
    }
}

void Knob::setToMax(){
    setRelativeValue(100.0);
}


double Knob::getRelativeValue() const{
    return _relativeValue;
}

double Knob::getRealValue() const{
    return _realValue;
}

bool Knob::autoFind() const{
    return getAllowedValues().size() > 1;
}

KnobWorkers::KnobWorkers(KnobConfWorkers confWorkers, ff_farm<>& farm):
        _confWorkers(confWorkers), _farm(farm){
    _emitter = dynamic_cast<AdaptiveNode*>(_farm.getEmitter());
    svector<ff_node*> workers = _farm.getNWorkers();
    for(size_t i = 0; i < workers.size(); i++){
        _activeWorkers.push_back(dynamic_cast<AdaptiveNode*>(workers[i]));
    }
    _collector = dynamic_cast<AdaptiveNode*>(_farm.getCollector());

    if(confWorkers == KNOB_WORKERS_YES){
        for(size_t i = 0; i < workers.size(); i++){
            _knobValues.push_back(i + 1);
        }
    }else{
        _knobValues.push_back(workers.size());
    }
}

void KnobWorkers::changeValueReal(double v){
    prepareToFreeze();
    freeze();

    changeActiveNodes(v);
    notifyNewConfiguration(v);

    prepareToRun(v);
    run(v);
}

std::vector<double> KnobWorkers::getAllowedValues() const{
    return _knobValues;
}

uint KnobWorkers::getNumActiveWorkers() const{
    return _activeWorkers.size();
}

const std::vector<AdaptiveNode*>& KnobWorkers::getActiveWorkers() const{
    return _activeWorkers;
}

const std::vector<AdaptiveNode*>& KnobWorkers::getInactiveWorkers() const{
    return _inactiveWorkers;
}

void KnobWorkers::prepareToFreeze(){
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

void KnobWorkers::freeze(){
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

void KnobWorkers::prepareToRun(uint numWorkers){
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

void KnobWorkers::run(uint numWorkers){
    if(_emitter){
        _emitter->thawAll(numWorkers);
    }
    if(_collector){
        _farm.getgt()->thaw(true, numWorkers);
    }
}

void KnobWorkers::changeActiveNodes(uint v){
    uint workersNumDiff = abs(_realValue - v);
    if(_realValue > v){
        /** Move workers from active to inactive. **/
        moveEndToFront(_activeWorkers, _inactiveWorkers, workersNumDiff);
    }else{
        /** Move workers from inactive to active. **/
        moveFrontToEnd(_inactiveWorkers, _activeWorkers, workersNumDiff);
    }
}

void KnobWorkers::notifyNewConfiguration(uint numWorkers){
    if(_emitter){
        _emitter->notifyWorkersChange(_realValue, numWorkers);
    }
    for(size_t i = 0; i < numWorkers; i++){
        AdaptiveNode* w = _activeWorkers.at(i);
        w->notifyWorkersChange(_realValue, numWorkers);
    }
    if(_collector){
        _collector->notifyWorkersChange(_realValue, numWorkers);
    }
}

KnobMapping::KnobMapping(KnobConfMapping confMapping,
                         KnobConfSNodeMapping confEmitterMapping,
                         KnobConfSNodeMapping confCollectorMapping,
                         KnobConfHyperthreading confHyperthreading,
                         const Mammut& mammut,
                         AdaptiveNode* emitter,
                         AdaptiveNode* collector,
                         const KnobWorkers& knobWorkers):
                             _confMapping(confMapping),
                             _confEmitterMapping(confEmitterMapping),
                             _confCollectorMapping(confCollectorMapping),
                             _confHyperthreading(confHyperthreading),
                             _emitterVirtualCore(NULL),
                             _collectorVirtualCore(NULL),
                             _emitter(emitter),
                             _collector(collector),
                             _knobWorkers(knobWorkers),
                             _topologyHandler(mammut.getInstanceTopology()){
    ;
}

void KnobMapping::changeValueReal(double v){
    setVcOrder();

    switch((KnobConfMapping) v){
        case KNOB_MAPPING_LINEAR:{
            performLinearMapping();
        }break;
        default:{
            throw runtime_error("KnobMapping: This situation should never happen.");
        }break;
    }

    /** Updates unused virtual cores. **/
    for(size_t i = 0; i < _vcOrder.size(); i++){
        VirtualCore* vc = _vcOrder.at(i);
        if(vc != _emitterVirtualCore && vc != _collectorVirtualCore &&
           !contains(_activeVirtualCores, vc) &&
           !contains(_inactiveVirtualCores, vc)){
            _unusedVirtualCores.push_back(vc);
        }
    }
}

std::vector<double> KnobMapping::getAllowedValues() const{
    std::vector<double> r;
    switch(_confMapping){
        case KNOB_MAPPING_NO:{
            ;
        }break;
        case KNOB_MAPPING_AUTO:{
            r.push_back((double) ((int) KNOB_MAPPING_LINEAR));
            r.push_back((double) ((int) KNOB_MAPPING_CACHE_EFFICIENT));
        }break;
        case KNOB_MAPPING_LINEAR:{
            r.push_back((double) ((int) KNOB_MAPPING_LINEAR));
        }break;
        case KNOB_MAPPING_CACHE_EFFICIENT:{
            r.push_back((double) ((int) KNOB_MAPPING_CACHE_EFFICIENT));
        }break;
    }

    return r;
}

mammut::topology::VirtualCore* KnobMapping::getEmitterVirtualCore() const{
    return _emitterVirtualCore;
}

mammut::topology::VirtualCore* KnobMapping::getCollectorVirtualCore() const{
    return _collectorVirtualCore;
}

const std::vector<mammut::topology::VirtualCore*>& KnobMapping::getWorkersVirtualCore() const{
    return _workersVirtualCores;
}

const std::vector<mammut::topology::VirtualCore*>& KnobMapping::getActiveVirtualCores() const{
    return _activeVirtualCores;
}

const std::vector<mammut::topology::VirtualCore*>& KnobMapping::getInactiveVirtualCores() const{
    return _inactiveVirtualCores;
}

const std::vector<mammut::topology::VirtualCore*>& KnobMapping::getUnusedVirtualCores() const{
    return _unusedVirtualCores;
}


void KnobMapping::getMappingIndexes(size_t& emitterIndex,
                                    size_t& firstWorkerIndex,
                                    size_t& collectorIndex){
    size_t nextIndex = 0;
    if(_emitter && !_emitterVirtualCore){
        emitterIndex = nextIndex;
        if(_confEmitterMapping == KNOB_SNODE_MAPPING_ALONE){
            nextIndex = (nextIndex + 1) % _vcOrder.size();
        }
    }

    firstWorkerIndex = nextIndex;
    nextIndex = (nextIndex + (uint)_knobWorkers.getRealValue()) %
                _vcOrder.size();
    //TODO Con troppi workers potrebbe non rispettare l'alone.
    // SOluzione: Invece di avere direttamente l'id, mettere gli id in un
    // vettore e rimuovere gli id gia utilizzati da emitter o collector
    if(_collector && !_collectorVirtualCore){
        if(_confCollectorMapping == KNOB_SNODE_MAPPING_COLLAPSED){
            nextIndex = (nextIndex - 1) % _vcOrder.size();
        }
        collectorIndex = nextIndex;
    }
}

void KnobMapping::setVcOrder(){
    switch(_confMapping){
        case KNOB_MAPPING_LINEAR:{
           /*
            * Generates a vector of virtual cores to be used for linear
            * mapping.node. It contains first one virtual core per physical
            * core (virtual cores on the same CPU are consecutive).
            * Then, the other groups of virtual cores follow.
            */
            vector<Cpu*> cpus = _topologyHandler->getCpus();

            size_t virtualUsed = 0;
            size_t virtualPerPhysical;
            if(_confHyperthreading != KNOB_HT_NO){
                virtualPerPhysical = _topologyHandler->getVirtualCores().size() /
                        _topologyHandler->getPhysicalCores().size();
            }else{
                virtualPerPhysical = 1;
            }
            while(virtualUsed < virtualPerPhysical){
                for(size_t i = 0; i < cpus.size(); i++){
                    vector<PhysicalCore*> phyCores = cpus.at(i)->
                                                     getPhysicalCores();
                    for(size_t j = 0; j < phyCores.size(); j++){
                        _vcOrder.push_back(phyCores.at(j)->getVirtualCores().
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
}

void KnobMapping::performLinearMapping(){
    size_t emitterIndex, firstWorkerIndex, collectorIndex;
    getMappingIndexes(emitterIndex, firstWorkerIndex, collectorIndex);

    const vector<AdaptiveNode*>& activeWorkers = _knobWorkers.getActiveWorkers();
    const vector<AdaptiveNode*>& inactiveWorkers = _knobWorkers.getInactiveWorkers();

    vector<VirtualCore*> remainingCores = _vcOrder;
    reverse(remainingCores.begin(), remainingCores.end());
    _activeVirtualCores.clear();
    _inactiveVirtualCores.clear();

    if(_emitter && _confEmitterMapping != KNOB_SNODE_MAPPING_NO){
        _emitterVirtualCore = _vcOrder.at(emitterIndex);
        _activeVirtualCores.push_back(_emitterVirtualCore);
        if(!_emitterVirtualCore->isHotPlugged()){
            _emitterVirtualCore->hotPlug();
        }
        _emitter->move(_emitterVirtualCore);
        remainingCores.pop_back();
    }


    for(size_t i = 0; i < activeWorkers.size(); i++){
        VirtualCore* vc = _vcOrder.at((firstWorkerIndex + i) %
                          _vcOrder.size());
        _workersVirtualCores.push_back(vc);
        _activeVirtualCores.push_back(vc);
        if(!vc->isHotPlugged()){
            vc->hotPlug();
        }
        activeWorkers.at(i)->move(vc);
        remainingCores.pop_back();
    }

    if(_collector && _confCollectorMapping != KNOB_SNODE_MAPPING_NO){
        _collectorVirtualCore = _vcOrder.at(collectorIndex);
        _activeVirtualCores.push_back(_collectorVirtualCore);
        if(!_collectorVirtualCore->isHotPlugged()){
            _collectorVirtualCore->hotPlug();
        }
        _collector->move(_collectorVirtualCore);
        remainingCores.pop_back();
    }

    for(size_t i = 0; i < inactiveWorkers.size() && !remainingCores.empty(); i++){
        _inactiveVirtualCores.push_back(remainingCores.back());
        remainingCores.pop_back();
    }
}

/*
void KnobMapping::updateUsedCpus(){
    _usedCpus.clear();
    _unusedCpus.clear();

    for(size_t i = 0; i < _activeVirtualCores.size(); i++){
        CpuId cpuId = _activeVirtualCores.at(i)->getCpuId();
        if(!contains(_usedCpus, cpuId)){
            _usedCpus.push_back(cpuId);
        }
    }

    vector<Cpu*> cpus = _topologyHandler->getCpus();
    for(size_t i = 0; i < cpus.size(); i++){
        if(!contains(_usedCpus, cpus.at(i)->getCpuId())){
            _unusedCpus.push_back(cpus.at(i)->getCpuId());
        }
    }
}
*/

KnobFrequency::KnobFrequency(KnobConfFrequencies confFrequency,
                             const Mammut& mammut, bool useTurboBoost,
                             StrategyUnusedVirtualCores inactiveVc,
                             StrategyUnusedVirtualCores unusedVc,
                             const KnobMapping& knobMapping):
        _confFrequency(confFrequency),
        _frequencyHandler(mammut.getInstanceCpuFreq()),
        _topologyHandler(mammut.getInstanceTopology()),
        _cpufreqHandle(mammut.getInstanceCpuFreq()),
        _inactiveVc(inactiveVc),
        _unusedVc(unusedVc),
        _knobMapping(knobMapping){
    _availableFrequencies = _frequencyHandler->getDomains().at(0)->getAvailableFrequencies();

    // Remove turbo boost frequency.
    if(!useTurboBoost && !_availableFrequencies.empty()){
        if(intToString(_availableFrequencies.back()).at(3) == '1'){
            _availableFrequencies.pop_back();
        }
    }

    if(_availableFrequencies.empty()){
        // Insert dummy constant frequency
        _availableFrequencies.push_back(1.0);
    }

    if(_confFrequency == KNOB_FREQUENCY_YES){
        for(size_t i = 0; i < _availableFrequencies.size(); i++){
            _allowedValues.push_back(_availableFrequencies.at(i));
        }
    }
}

void KnobFrequency::changeValueReal(double v){
    _scalableDomains = _frequencyHandler->getDomains(_knobMapping.getActiveVirtualCores());
    Domain* currentDomain;
    for(size_t i = 0; i < _scalableDomains.size(); i++){
        currentDomain = _scalableDomains.at(i);
        if(!currentDomain->setFrequencyUserspace((uint)v)){
            throw runtime_error("AdaptivityManagerFarm: Impossible "
                                "to set the specified frequency.");
        }
    }
}

std::vector<double> KnobFrequency::getAllowedValues() const{
    return _allowedValues;
}

void KnobFrequency::applyUnusedVCStrategyOff(const vector<VirtualCore*>& unusedVc){
    for(size_t i = 0; i < unusedVc.size(); i++){
        VirtualCore* vc = unusedVc.at(i);
        if(vc->isHotPluggable() && vc->isHotPlugged()){
            vc->hotUnplug();
        }
    }
}

void KnobFrequency::applyUnusedVCStrategyLowestFreq(const vector<VirtualCore*>& vc){
    vector<Domain*> unusedDomains = _cpufreqHandle->getDomainsComplete(vc);
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

void KnobFrequency::applyUnusedVCStrategy(){
    /**
     * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the
     * virtual cores on a domain, we can also lower its frequency to
     * the minimum.
     */
    vector<VirtualCore*> inactiveVc = _knobMapping.getInactiveVirtualCores();
    vector<VirtualCore*> unusedVc = _knobMapping.getUnusedVirtualCores();

    vector<VirtualCore*> virtualCores;
    if(_inactiveVc != STRATEGY_UNUSED_VC_NONE){
        insertToEnd(inactiveVc, virtualCores);
    }
    if(_unusedVc != STRATEGY_UNUSED_VC_NONE){
        insertToEnd(unusedVc, virtualCores);
    }
    applyUnusedVCStrategyLowestFreq(virtualCores);


    virtualCores.clear();
    if(_inactiveVc == STRATEGY_UNUSED_VC_OFF){
        insertToEnd(inactiveVc, virtualCores);
    }
    if(_unusedVc == STRATEGY_UNUSED_VC_OFF){
        insertToEnd(unusedVc, virtualCores);
    }
    applyUnusedVCStrategyOff(virtualCores);
}

}


