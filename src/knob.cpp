/*
 * knob.hpp
 *
 * Created on: 02/11/2015
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
 */

#include "knob.hpp"
#include "parameters.hpp"

#include "external/Mammut/mammut/mammut.hpp"
#include "external/Mammut/mammut/utils.hpp"
#include "external/Mammut/mammut/cpufreq/cpufreq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_KNOB
#define DEBUG(x) do { cerr << "[Knob] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

using namespace ff;

using namespace std;

std::string knobTypeToString(KnobType kv){
    switch(kv){
        case KNOB_TYPE_VIRTUAL_CORES:{
            return "Cores";
        }break;
        case KNOB_TYPE_HYPERTHREADING:{
            return "HTLevel";
        }break;
        case KNOB_TYPE_MAPPING:{
            return "Mapping";
        }break;
        case KNOB_TYPE_FREQUENCY:{
            return "Frequency";
        }break;
        default:{
            return "Unknown";
        }break;
    }
}

bool Knob::getRealFromRelative(double relative, double& real) const{
    // Maps from the range [0, 100] to the real range.
    vector<double> values = getAllowedValues();
    if(values.size()){
        uint index = round((double)(values.size() - 1) * (relative / 100.0));
        real = values.at(index);
        return true;
    }else{
        return false;
    }
}

void Knob::setRelativeValue(double v){
    double real;
    if(getRealFromRelative(v, real)){
        setRealValue(real);
    }
}

void Knob::setRealValue(double v){
    if(getAllowedValues().size()){
        changeValueReal(v);
        _realValue = v;
    }
}

void Knob::setToMax(){
    setRelativeValue(100.0);
}

void Knob::lock(double v){
    setRealValue(v);
    _knobValues.clear();
    _knobValues.push_back(v);
    _locked = true;
}

void Knob::lockToMax(){
    if(_knobValues.size()){
        lock(*_knobValues.end());
    }else{
        _locked = true;
    }
}

void Knob::lockToMin(){
    if(_knobValues.size()){
        lock(*_knobValues.begin());
    }else{
        _locked = true;
    }
}

bool Knob::isLocked() const{
    return _locked;
}

double Knob::getRealValue() const{
    return _realValue;
}

std::vector<double> Knob::getAllowedValues() const{
    return _knobValues;
}

KnobVirtualCores::KnobVirtualCores(Parameters p, uint plus):_p(p), _plus(plus){
    std::vector<VirtualCore*> virtualCores = _p.mammut.getInstanceTopology()->getVirtualCores();
    changeMax(virtualCores.size() - _plus);
}

void KnobVirtualCores::changeValueReal(double v){;}

void KnobVirtualCores::changeMax(double v){
    _knobValues.clear();
    for(size_t i = 0; i < v; i++){
        if(!utils::contains(_p.disallowedNumCores, (uint) i + 1)){
            _knobValues.push_back(i + 1);
        }
    }
}

#ifdef DEBUG_KNOB
std::ostream& operator<< (std::ostream& out, const std::vector<AdaptiveNode*>& v){
    out << "[";
    for(size_t i = 0; i < v.size(); i++){
        out << (v.at(i))->getOSThreadId() << ", ";
    }
    out << "]";
    return out;
}

std::ostream& operator<< (std::ostream& out, const std::vector<VirtualCore*>& v){
    out << "[";
    for(size_t i = 0; i < v.size(); i++){
        out << (v.at(i))->getVirtualCoreId() << ", ";
    }
    out << "]";
    return out;
}

std::ostream& operator<< (std::ostream& out, const std::vector<Domain*>& v){
    out << "[";
    for(size_t i = 0; i < v.size(); i++){
        out << (v.at(i))->getId() << ", ";
    }
    out << "]";
    return out;
}
#endif

KnobVirtualCoresFarm::KnobVirtualCoresFarm(Parameters p,
                             AdaptiveNode* emitter, AdaptiveNode* collector,
                             ff::ff_gatherer* gt,
                             const std::vector<AdaptiveNode*> workers,
                             const volatile bool* terminated):
                KnobVirtualCores(p, (emitter?1:0) + (collector?1:0)),
        _emitter(emitter), _collector(collector),
        _gt(gt), _allWorkers(workers), _terminated(terminated){
    _realValue = _allWorkers.size();
    _knobValues.clear();
    for(size_t i = 0; i < _allWorkers.size(); i++){
        if(!utils::contains(p.disallowedNumCores, (uint) i + 1)){
            _knobValues.push_back(i + 1);
        }
    }

    _activeWorkers = _allWorkers;
    DEBUG("Knob workers created.");
}

void KnobVirtualCoresFarm::changeValueReal(double v){
    if(v != _realValue){
        DEBUG("[Workers] Changing real value to: " << v);
        prepareToFreeze();
        freeze();

        if(!*_terminated){
            _activeWorkers = vector<AdaptiveNode*>(_allWorkers.begin(),
                                                   _allWorkers.begin() + v);

            notifyNewConfiguration(v);

            prepareToRun(v);
            run(v);
            DEBUG("[Workers] Active Workers: " << _activeWorkers);
        }
    }
}

std::vector<double> KnobVirtualCoresFarm::getAllowedValues() const{
    return _knobValues;
}

std::vector<AdaptiveNode*> KnobVirtualCoresFarm::getActiveWorkers() const{
    return _activeWorkers;
}

void KnobVirtualCoresFarm::prepareToFreeze(){
    DEBUG("[Workers] Preparing the farm for freezing.");
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

void KnobVirtualCoresFarm::freeze(){
    DEBUG("[Workers] Freezing the farm.");
    if(_emitter){
        _emitter->freezeAll((void*)(_collector?FF_EOSW:FF_GO_OUT));
    }

    for(size_t i = 0; i < _activeWorkers.size(); i++){
        DEBUG("Waiting Worker " << i << " to be frozen:");
        _activeWorkers.at(i)->wait_freezing();
        assert(_activeWorkers.at(i)->isfrozen());
        DEBUG("Worker " << i << " frozen:");
    }
    if(_collector){
        _gt->wait_freezing();
    }
}

void KnobVirtualCoresFarm::prepareToRun(uint numWorkers){
    DEBUG("[Workers] Running the farm.");
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

void KnobVirtualCoresFarm::run(uint numWorkers){
    DEBUG("[Workers] Running the emitter.");
    if(_emitter){
        _emitter->thawAll(numWorkers);
    }
    DEBUG("[Workers] Running the collector.");
    if(_collector){
        _gt->thaw(true, numWorkers);
    }
}

void KnobVirtualCoresFarm::notifyNewConfiguration(uint numWorkers){
    if(_emitter){
        _emitter->notifyRethreading(_realValue, numWorkers);
    }
    for(size_t i = 0; i < numWorkers; i++){
        AdaptiveNode* w = _activeWorkers.at(i);
        w->notifyRethreading(_realValue, numWorkers);
    }
    if(_collector){
        _collector->notifyRethreading(_realValue, numWorkers);
    }
}


KnobHyperThreading::KnobHyperThreading(Parameters p){
    vector<PhysicalCore*> physical = p.mammut.getInstanceTopology()->getPhysicalCores();
    size_t maxHtLevel = physical.at(0)->getVirtualCores().size();
    for(size_t i = 0; i < maxHtLevel; i++){
        _knobValues.push_back(i + 1);
    }
}

void KnobHyperThreading::changeValueReal(double v){;}

KnobMapping::KnobMapping(const Parameters& p,
                         const KnobVirtualCores& knobCores,
                         const KnobHyperThreading& knobHyperThreading):
         _p(p),
         _knobCores(knobCores),
         _knobHyperThreading(knobHyperThreading),
         _topologyHandler(p.mammut.getInstanceTopology()){
    for(size_t i = 0; i < MAPPING_TYPE_NUM; i++){
        _knobValues.push_back((MappingType) i);
    }
    _realValue = MAPPING_TYPE_LINEAR; // This is just for initialization, is not locked.
}

template<> char const* enumStrings<MappingType>::data[] = {
    "LINEAR",
    "INTERLEAVED",
    "CACHE_OPTIMAL"
};

void KnobMapping::changeValueReal(double v){
    DEBUG("[Mapping] Changing real value to: " << enumToString<MappingType>((MappingType)v));

    // Its length will be equal to _knobCores.getRealValue()
    vector<VirtualCore*> vcOrder;
    switch((MappingType) v){
        case MAPPING_TYPE_LINEAR:{
            vcOrder = computeVcOrderLinear();
        }break;
        case MAPPING_TYPE_INTERLEAVED:{
            vcOrder = computeVcOrderInterleaved();
        }break;
        default:{
            throw runtime_error("KnobMapping: Mapping type still not supported.");
        }break;
    }

    /** Performs mapping. **/
    _activeVirtualCores = vcOrder;
    move(vcOrder);

    /** Updates unused virtual cores. **/
    _unusedVirtualCores.clear();
    vector<VirtualCore*> allVcs = _topologyHandler->getVirtualCores();
    for(size_t i = 0; i < allVcs.size(); i++){
        VirtualCore* vc = allVcs.at(i);
        if(!contains(_activeVirtualCores, vc)){
            _unusedVirtualCores.push_back(vc);
        }
    }

    DEBUG("[Mapping] Active VCs: " << _activeVirtualCores);
    DEBUG("[Mapping] Unused VCs: " << _unusedVirtualCores);
}

void KnobMapping::setAllowedCores(std::vector<mammut::topology::VirtualCore*> vc){
    _allowedVirtualCores = vc;
}

bool KnobMapping::isAllowed(mammut::topology::VirtualCore* v) const{
    if(!_allowedVirtualCores.size()){
        return true;
    }else{
        return utils::contains(_allowedVirtualCores, v);
    }
}

const vector<VirtualCore*>& KnobMapping::getActiveVirtualCores() const{
    return _activeVirtualCores;
}

const vector<VirtualCore*>& KnobMapping::getUnusedVirtualCores() const{
    return _unusedVirtualCores;
}

size_t KnobMapping::getNumVirtualCores(){
    return _knobCores.getRealValue();
}

vector<VirtualCore*> KnobMapping::computeVcOrderLinear(){
   /*
    * Generates a vector of virtual cores to be used for linear
    * mapping. It contains first one virtual core per physical
    * core (virtual cores on the same CPU are consecutive).
    * Then, the other groups of virtual cores follow.
    */
    vector<VirtualCore*> vcOrder;
    size_t virtualPerPhysical = _knobHyperThreading.getRealValue();

    vector<Cpu*> cpus = _topologyHandler->getCpus();
    for(size_t i = 0; i < cpus.size(); i++){
        vector<PhysicalCore*> phyCores = cpus.at(i)->getPhysicalCores();
        for(size_t j = 0; j < phyCores.size(); j++){
            vector<VirtualCore*> virtCores = phyCores.at(j)->getVirtualCores();
            for(size_t k = 0; k < virtCores.size(); k++){
                if(k >= virtualPerPhysical){break;}
                if(vcOrder.size() == getNumVirtualCores()){
                    return vcOrder;
                }
                if((!_p.isolateManager || virtCores.at(k)->getVirtualCoreId() != MANAGER_VIRTUAL_CORE) &&
                    isAllowed(virtCores.at(k))){
                    vcOrder.push_back(virtCores.at(k));
                }
            }
        }
    }
    return vcOrder;
}

vector<VirtualCore*> KnobMapping::computeVcOrderInterleaved(){
   /*
    * Generates a vector of virtual cores to be used for interleaved
    * mapping.
    */
    vector<VirtualCore*> vcOrder;
    size_t virtualPerPhysical = _knobHyperThreading.getRealValue();

    vector<Cpu*> cpus = _topologyHandler->getCpus();
    size_t nextPhysicalId = 0;
    size_t nextVirtualId = 0;
    while(true){
        for(size_t i = 0; i < cpus.size(); i++){
            vector<PhysicalCore*> phyCores = cpus.at(i)->getPhysicalCores();
            PhysicalCore* p = phyCores.at(nextPhysicalId);
            vector<VirtualCore*> virtCores = p->getVirtualCores();
            VirtualCore* v = virtCores.at(nextVirtualId);
            if((!_p.isolateManager || v->getVirtualCoreId() != MANAGER_VIRTUAL_CORE) &&
                isAllowed(v)){
                vcOrder.push_back(v);
                if(vcOrder.size() == getNumVirtualCores()){return vcOrder;}
            }
        }
        if(nextPhysicalId + 1 < cpus.at(0)->getPhysicalCores().size()){
            nextPhysicalId += 1;
        }else{
            nextPhysicalId = 0;
            nextVirtualId += 1;
            if(nextVirtualId >= virtualPerPhysical){return vcOrder;}
        }
    }
}

KnobMappingExternal::KnobMappingExternal(const Parameters& p,
            const KnobVirtualCores& knobCores,
            const KnobHyperThreading& knobHyperThreading):
        KnobMapping(p, knobCores, knobHyperThreading), _processHandler(NULL){
    ;
}

void KnobMappingExternal::setPid(pid_t pid){
    setProcessHandler(_p.mammut.getInstanceTask()->getProcessHandler(pid));
}

void KnobMappingExternal::setProcessHandler(task::ProcessHandler* processHandler){
    _processHandler = processHandler;
}

void KnobMappingExternal::move(const vector<VirtualCore*>& vcOrder){
    if(_processHandler){
        _processHandler->move(vcOrder);
    }
}

KnobMappingFarm::KnobMappingFarm(const Parameters& p,
            const KnobVirtualCoresFarm& knobCores,
            const KnobHyperThreading& knobHyperThreading,
            AdaptiveNode* emitter,
            AdaptiveNode* collector):
        KnobMapping(p, knobCores, knobHyperThreading), _emitter(emitter),
        _collector(collector){
    ;
}

size_t KnobMappingFarm::getNumVirtualCores(){
    size_t v = _knobCores.getRealValue();
    if(_emitter) ++v;
    if(_collector) ++v;
    return v;
}

void KnobMappingFarm::move(const vector<VirtualCore*>& vcOrder){
    vector<AdaptiveNode*> workers = ((KnobVirtualCoresFarm*) &_knobCores)->getActiveWorkers();
    size_t nextIndex = 0;
    if(_emitter){
        _emitter->move((VirtualCore*) vcOrder[nextIndex]);
        nextIndex = (nextIndex + 1) % vcOrder.size();
    }

    for(size_t i = 0; i < workers.size(); i++){
        workers[i]->move((VirtualCore*) vcOrder[nextIndex]);
        nextIndex = (nextIndex + 1) % vcOrder.size();
    }

    if(_collector){
        _collector->move((VirtualCore*) vcOrder[nextIndex]);
        nextIndex = (nextIndex + 1) % vcOrder.size();
    }
}


KnobFrequency::KnobFrequency(Parameters p, const KnobMapping& knobMapping):
        _p(p),
        _knobMapping(knobMapping),
        _frequencyHandler(_p.mammut.getInstanceCpuFreq()),
        _topologyHandler(_p.mammut.getInstanceTopology()),
        _cpufreqHandle(_p.mammut.getInstanceCpuFreq()){

    std::vector<mammut::cpufreq::Frequency> availableFrequencies;
    availableFrequencies = _frequencyHandler->getDomains().at(0)->getAvailableFrequencies();
    _realValue = availableFrequencies.front();

    std::vector<mammut::cpufreq::Domain*> scalableDomains;
    scalableDomains = _frequencyHandler->getDomains();
    Domain* currentDomain;
    for(size_t i = 0; i < scalableDomains.size(); i++){
        currentDomain = scalableDomains.at(i);
        if(!currentDomain->setGovernor(GOVERNOR_USERSPACE)){
            throw runtime_error("AdaptivityManagerFarm: Impossible "
                                "to set the specified governor.");
        }
    }
    for(size_t i = 0; i < availableFrequencies.size(); i++){
        _knobValues.push_back(availableFrequencies.at(i));
    }
}

void KnobFrequency::changeValueReal(double v){
    DEBUG("[Frequency] Changing real value to: " << v);
    std::vector<mammut::cpufreq::Domain*> scalableDomains;
    scalableDomains = _frequencyHandler->getDomains(_knobMapping.getActiveVirtualCores());
    Domain* currentDomain;
    for(size_t i = 0; i < scalableDomains.size(); i++){
        currentDomain = scalableDomains.at(i);
        if(!currentDomain->setFrequencyUserspace((uint)v)){
            throw runtime_error("AdaptivityManagerFarm: Impossible "
                                "to set the specified frequency.");
        }
    }
    applyUnusedVCStrategy(v);
    DEBUG("[Frequency] Frequency changed for domains: " << scalableDomains);
    DEBUG("[Frequency] Active VC: " << _knobMapping.getActiveVirtualCores());
    DEBUG("[Frequency] Unused VC: " << _knobMapping.getUnusedVirtualCores());
}

void KnobFrequency::applyUnusedVCStrategySame(const vector<VirtualCore*>& unusedVc, Frequency v){
    vector<Domain*> unusedDomains = _cpufreqHandle->getDomainsComplete(unusedVc);
    DEBUG("[Frequency] " << unusedDomains.size() << " unused domains.");
    for(size_t i = 0; i < unusedDomains.size(); i++){
        Domain* domain = unusedDomains.at(i);
        DEBUG("[Frequency] Setting unused domain " << domain->getId() << " to: " << v);
        if(!domain->setFrequencyUserspace((uint)v)){
            throw runtime_error("AdaptivityManagerFarm: Impossible "
                                "to set the specified frequency.");
        }    
    }
}


void KnobFrequency::applyUnusedVCStrategyOff(const vector<VirtualCore*>& unusedVc){
    for(size_t i = 0; i < unusedVc.size(); i++){
        VirtualCore* vc = unusedVc.at(i);
        if(vc->isHotPluggable() && vc->isHotPlugged()){
            vc->hotUnplug();
        }
    }
}

void KnobFrequency::applyUnusedVCStrategyLowestFreq(const vector<VirtualCore*>& unusedVc){
    vector<Domain*> unusedDomains = _cpufreqHandle->getDomainsComplete(unusedVc);
    for(size_t i = 0; i < unusedDomains.size(); i++){
        Domain* domain = unusedDomains.at(i);
        if(!domain->setGovernor(GOVERNOR_USERSPACE) ||
           !domain->setLowestFrequencyUserspace()){
            throw runtime_error("AdaptivityManagerFarm: Impossible to "
                                "set lowest frequency for unused "
                                "virtual cores.");
        }
    }
}

void KnobFrequency::applyUnusedVCStrategy(Frequency v){
    /**
     * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the
     * virtual cores on a domain, we can also lower its frequency to
     * the minimum.
     */
    vector<VirtualCore*> unusedVc = _knobMapping.getUnusedVirtualCores();

    if(_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_SAME){
        applyUnusedVCStrategySame(unusedVc, v);
    }else{
        vector<VirtualCore*> virtualCores;
        if(_p.strategyUnusedVirtualCores != STRATEGY_UNUSED_VC_NONE){
            insertToEnd(unusedVc, virtualCores);
        }
        applyUnusedVCStrategyLowestFreq(virtualCores);

        virtualCores.clear();
        if(_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_OFF){
            insertToEnd(unusedVc, virtualCores);
        }
        applyUnusedVCStrategyOff(virtualCores);
    }
}

}


