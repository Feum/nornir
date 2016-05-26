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
        case KNOB_TYPE_WORKERS:{
            return "Workers";
        }break;
        case KNOB_TYPE_FREQUENCY:{
            return "Frequency";
        }break;
        case KNOB_TYPE_MAPPING:{
            return "Mapping";
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

double Knob::getRealValue() const{
    return _realValue;
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
#endif

KnobWorkers::KnobWorkers(Parameters p, AdaptiveNode* emitter,
                         AdaptiveNode* collector, ff::ff_gatherer* gt,
                         const std::vector<AdaptiveNode*> workers,
                         const volatile bool* terminated):
        _confWorkers(p.knobWorkers), _emitter(emitter), _collector(collector),
        _gt(gt), _allWorkers(workers), _terminated(terminated){
    _realValue = _allWorkers.size();

    if(needsCalibration()){
        for(size_t i = 0; i < _allWorkers.size(); i++){
            if(!utils::contains(p.disallowedNumCores, (uint) i + 1)){
                _knobValues.push_back(i + 1);
            }
        }
    }else{
        _knobValues.push_back(_allWorkers.size());
    }
    _activeWorkers = _allWorkers;
    DEBUG("Knob workers created.");
}

bool KnobWorkers::needsCalibration() const{
    return _confWorkers != KNOB_WORKERS_NO;
}

void KnobWorkers::changeValueReal(double v){
    if(v != _realValue && _confWorkers == KNOB_WORKERS_THREADS){
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

std::vector<double> KnobWorkers::getAllowedValues() const{
    return _knobValues;
}

uint KnobWorkers::getWorkersPhysicalCores() const{
    switch(_confWorkers){
        case KNOB_WORKERS_NO:
        case KNOB_WORKERS_THREADS:{
            return _activeWorkers.size();
        }break;
        case KNOB_WORKERS_MAPPING:{
            // _activeWorkers.size() is the number of active threads.
            // However, it never changes since we map more threads
            // to fewer cores, so we return the real value of the knob.
            return getRealValue();
        }break;
        default:{
            throw runtime_error("Unkown KnobWorkers configuration parameter.");
        }
    }
}

const std::vector<AdaptiveNode*>& KnobWorkers::getActiveWorkers() const{
    return _activeWorkers;
}

void KnobWorkers::prepareToFreeze(){
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

void KnobWorkers::freeze(){
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

void KnobWorkers::prepareToRun(uint numWorkers){
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

void KnobWorkers::run(uint numWorkers){
    DEBUG("[Workers] Running the emitter.");
    if(_emitter){
        _emitter->thawAll(numWorkers);
    }
    DEBUG("[Workers] Running the collector.");
    if(_collector){
        _gt->thaw(true, numWorkers);
    }
}

void KnobWorkers::notifyNewConfiguration(uint numWorkers){
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

KnobMapping::KnobMapping(const Parameters& p,
                         AdaptiveNode* emitter,
                         AdaptiveNode* collector,
                         const KnobWorkers& knobWorkers):
                             _p(p),
                             _vcOrder(_vcOrderCacheEfficient),
                             _emitterVirtualCore(NULL),
                             _collectorVirtualCore(NULL),
                             _emitter(emitter),
                             _collector(collector),
                             _knobWorkers(knobWorkers),
                             _topologyHandler(p.mammut.getInstanceTopology()){
    computeVcOrderLinear();
    _realValue = KNOB_MAPPING_LINEAR;
}

bool KnobMapping::needsCalibration() const{
    return _p.knobMapping == KNOB_MAPPING_AUTO;
}

#ifdef DEBUG_KNOB
std::ostream& operator<< (std::ostream& out, const std::vector<VirtualCore*>& v){
    out << "[";
    for(size_t i = 0; i < v.size(); i++){
        out << (v.at(i))->getVirtualCoreId() << ", ";
    }
    out << "]";
    return out;
}
#endif

void KnobMapping::changeValueReal(double v){
    DEBUG("[Mapping] Changing real value to: " << enumToString<KnobConfMapping>((KnobConfMapping)v));

    switch((KnobConfMapping) v){
        case KNOB_MAPPING_LINEAR:{
            _vcOrder = _vcOrderLinear;
            performLinearMapping();
        }break;
        default:{
            throw runtime_error("KnobMapping: This situation should never happen.");
        }break;
    }

    /** Updates unused virtual cores. **/
    _unusedVirtualCores.clear();
    vector<VirtualCore*> allVcs = _topologyHandler->getVirtualCores();
    for(size_t i = 0; i < allVcs.size(); i++){
        VirtualCore* vc = allVcs.at(i);
        if(vc != _emitterVirtualCore && vc != _collectorVirtualCore &&
           !contains(_activeVirtualCores, vc)){
            _unusedVirtualCores.push_back(vc);
        }
    }

    DEBUG("[Mapping] Active VCs: " << _activeVirtualCores);
    DEBUG("[Mapping] Unused VCs: " << _unusedVirtualCores);
}

std::vector<double> KnobMapping::getAllowedValues() const{
    std::vector<double> r;
    switch(_p.knobMapping){
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

const std::vector<mammut::topology::VirtualCore*>& KnobMapping::getUnusedVirtualCores() const{
    return _unusedVirtualCores;
}

void KnobMapping::computeVcOrderLinear(){
   /*
    * Generates a vector of virtual cores to be used for linear
    * mapping. It contains first one virtual core per physical
    * core (virtual cores on the same CPU are consecutive).
    * Then, the other groups of virtual cores follow.
    */
    vector<Cpu*> cpus = _topologyHandler->getCpus();

    size_t virtualUsed = 0;
    size_t virtualPerPhysical;
    if(_p.knobHyperthreading != KNOB_HT_NO){
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
                if(!_p.isolateManager || !phyCores.at(j)->getVirtualCore(MANAGER_VIRTUAL_CORE)){
                    _vcOrderLinear.push_back(phyCores.at(j)->getVirtualCores().
                                               at(virtualUsed));
                }
            }
        }
        ++virtualUsed;
    }
}

void KnobMapping::getMappingIndexes(size_t& emitterIndex,
                                    size_t& firstWorkerIndex,
                                    size_t& collectorIndex){
    size_t nextIndex = 0;
    if(_emitter){
        emitterIndex = nextIndex;
        if(_p.knobMappingEmitter == KNOB_SNODE_MAPPING_ALONE){
            nextIndex = (nextIndex + 1) % _vcOrder.size();
        }
    }

    if(_collector){
        collectorIndex = nextIndex;
        if(_p.knobMappingCollector == KNOB_SNODE_MAPPING_ALONE){
            nextIndex = (nextIndex + 1) % _vcOrder.size();
        }
    }

    firstWorkerIndex = nextIndex;
}

void KnobMapping::performLinearMapping(){
    size_t emitterIndex, firstWorkerIndex, collectorIndex;
    size_t numPhysicalCores = _topologyHandler->getPhysicalCores().size();
    getMappingIndexes(emitterIndex, firstWorkerIndex, collectorIndex);

    const vector<AdaptiveNode*>& activeWorkers = _knobWorkers.getActiveWorkers();

    _activeVirtualCores.clear();
    _workersVirtualCores.clear();

    /*********** Emitter mapping. ***********/
    if(_emitter && _p.knobMappingEmitter != KNOB_SNODE_MAPPING_NO){
        _emitterVirtualCore = _vcOrderLinear.at(emitterIndex);
        _activeVirtualCores.push_back(_emitterVirtualCore);
        if(!_emitterVirtualCore->isHotPlugged()){
            _emitterVirtualCore->hotPlug();
        }
        _emitter->move(_emitterVirtualCore);
    }

    /*********** Workers mapping. ***********/
    switch(_p.knobWorkers){
        case KNOB_WORKERS_MAPPING:{
            // Move workers over a set of cores.
            vector<const VirtualCore*> workersCores;
            for(size_t i = 0; i < _vcOrderLinear.size(); i++){
                if((!_emitter || i != emitterIndex) &&
                   (!_collector || i != collectorIndex)){
                    workersCores.push_back(_vcOrderLinear.at(i));
                }
            }

            for(size_t i = 0; i < activeWorkers.size(); i++){
                activeWorkers.at(i)->move(workersCores);
            }
        }break;
        case KNOB_WORKERS_THREADS:{
            size_t nextWorkerIndex = firstWorkerIndex;
            for(size_t i = 0; i < activeWorkers.size(); i++){
                VirtualCore* vc = _vcOrderLinear.at(nextWorkerIndex);

                _workersVirtualCores.push_back(vc);
                _activeVirtualCores.push_back(vc);
                if(!vc->isHotPlugged()){
                    vc->hotPlug();
                }

                activeWorkers.at(i)->move(vc);

                if(++nextWorkerIndex == _vcOrder.size()){
                    nextWorkerIndex = firstWorkerIndex;
                }
            }
        }break;
        default:{
            throw std::runtime_error("Unknown workers knob configuration.");
        }break;
    }

    /*********** Collector mapping. ***********/
    if(_collector && _p.knobMappingCollector != KNOB_SNODE_MAPPING_NO){
        _collectorVirtualCore = _vcOrderLinear.at(collectorIndex);
        _activeVirtualCores.push_back(_collectorVirtualCore);
        if(!_collectorVirtualCore->isHotPlugged()){
            _collectorVirtualCore->hotPlug();
        }
        _collector->move(_collectorVirtualCore);
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

#ifdef DEBUG_KNOB
std::ostream& operator<< (std::ostream& out, const std::vector<Domain*>& v){
    out << "[";
    for(size_t i = 0; i < v.size(); i++){
        out << (v.at(i))->getId() << ", ";
    }
    out << "]";
    return out;
}
#endif

KnobFrequency::KnobFrequency(KnobConfFrequencies confFrequency,
                             const Mammut& mammut, bool useTurboBoost,
                             StrategyUnusedVirtualCores unusedVc,
                             const KnobMapping& knobMapping):
        _confFrequency(confFrequency),
        _frequencyHandler(mammut.getInstanceCpuFreq()),
        _topologyHandler(mammut.getInstanceTopology()),
        _cpufreqHandle(mammut.getInstanceCpuFreq()),
        _unusedVc(unusedVc),
        _knobMapping(knobMapping){

    if(_confFrequency == KNOB_FREQUENCY_YES){
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
            _allowedValues.push_back(availableFrequencies.at(i));

        }
        DEBUG("[Frequency] Setting userspace governor. Scalable domains: " << scalableDomains);
        DEBUG("[Frequency] Unused VC strategy: " << _unusedVc);
    }
}

bool KnobFrequency::needsCalibration() const{
    return _confFrequency == KNOB_FREQUENCY_YES;
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
}

std::vector<double> KnobFrequency::getAllowedValues() const{
    return _allowedValues;
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

void KnobFrequency::applyUnusedVCStrategy(Frequency v){
    /**
     * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the
     * virtual cores on a domain, we can also lower its frequency to
     * the minimum.
     */
    vector<VirtualCore*> unusedVc = _knobMapping.getUnusedVirtualCores();

    if(_unusedVc == STRATEGY_UNUSED_VC_SAME){
        applyUnusedVCStrategySame(unusedVc, v);
    }else{
        vector<VirtualCore*> virtualCores;
        if(_unusedVc != STRATEGY_UNUSED_VC_NONE){
            insertToEnd(unusedVc, virtualCores);
        }
        applyUnusedVCStrategyLowestFreq(virtualCores);

        virtualCores.clear();
        if(_unusedVc == STRATEGY_UNUSED_VC_OFF){
            insertToEnd(unusedVc, virtualCores);
        }
        applyUnusedVCStrategyOff(virtualCores);
    }
}

}


