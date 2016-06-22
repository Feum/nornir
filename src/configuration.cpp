/*
 * configuration.cpp
 *
 * Created on: 05/12/2015
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

#include "configuration.hpp"

#ifdef DEBUG_CONFIGURATION
#define DEBUG(x) do { cerr << "[Configuration] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

Configuration::Configuration(const Parameters& p):
    _p(p), _combinationsCreated(false), _knobsChangeNeeded(false){
    ;
}

Configuration::~Configuration(){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        delete _knobs[i];
    }
}

//TODO: Works even if a std::vector is empty? (i.e. a knob has no values)
void Configuration::combinations(std::vector<std::vector<double> > array, size_t i, std::vector<double> accum){
    if(i == array.size()){
        KnobsValues kv(KNOB_VALUE_REAL);
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            kv[(KnobType) i] = accum.at(i);
        }
        _combinations.push_back(kv);
    }else{
        std::vector<double> row = array.at(i);
        for(size_t j = 0; j < row.size(); ++j){
            std::vector<double> tmp(accum);
            tmp.push_back(row[j]);
            combinations(array, i+1, tmp);
        }
    }
}

bool Configuration::equal(KnobsValues values) const{
    if(values.areReal()){
        return values == getRealValues();
    }else{
        double real;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            if(_knobs[(KnobType) i]->getRealFromRelative(values[(KnobType) i], real) &&
               real != getRealValue((KnobType) i)){
                return false;
            }
        }
        return true;
    }
}

bool Configuration::knobsChangeNeeded() const{
    return _knobsChangeNeeded;
}

void Configuration::createAllRealCombinations(){
    std::vector<std::vector<double>> values;
    std::vector<double> accum;
    _combinations.clear();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        values.push_back(_knobs[i]->getAllowedValues());
        if(!_knobs[i]->isLocked()){
            _knobsChangeNeeded = true;
        }
    }
    combinations(values, 0, accum);
    _combinationsCreated = true;
}

const std::vector<KnobsValues>& Configuration::getAllRealCombinations() const{
    if(!_combinationsCreated){
        throw std::runtime_error("[configuration.cpp] Combinations not created yet.");
    }
    return _combinations;
}

void Configuration::setFastReconfiguration(){
    if(_p.fastReconfiguration){
        ((KnobFrequency*) _knobs[KNOB_TYPE_FREQUENCY])->setRelativeValue(100.0);
    }
}

Knob* Configuration::getKnob(KnobType t) const{
    return _knobs[t];
}

void Configuration::maxAllKnobs(){
    if(_p.contractType == CONTRACT_NONE || !_knobsChangeNeeded){
        return;
    }
    DEBUG("Maxing all the knobs.");
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[(KnobType) i]->setToMax();
    }
}

double Configuration::getRealValue(KnobType t) const{
    return _knobs[t]->getRealValue();
}

KnobsValues Configuration::getRealValues() const{
    KnobsValues kv(KNOB_VALUE_REAL);
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        kv[(KnobType) i] = getRealValue((KnobType) i);
    }
    return kv;
}

bool Configuration::virtualCoresWillChange(const KnobsValues& values) const{
    double newNumWorkers = 0;
    if(values.areRelative()){
        assert(_knobs[KNOB_TYPE_VIRTUAL_CORES]->getRealFromRelative(values[KNOB_TYPE_VIRTUAL_CORES], newNumWorkers));
    }else{
        newNumWorkers = values[KNOB_TYPE_VIRTUAL_CORES];
    }
    return newNumWorkers != _knobs[KNOB_TYPE_VIRTUAL_CORES]->getRealValue();
}

ticks Configuration::startReconfigurationStatsKnob() const{
    if(_p.statsReconfiguration){
        return getticks();
    }
    return 0;
}

ticks Configuration::startReconfigurationStatsTotal() const{
    return startReconfigurationStatsKnob();
}

void Configuration::stopReconfigurationStatsKnob(ticks start, KnobType type,
                                                  bool vcChanged){
    if(_p.statsReconfiguration){
        if(type == KNOB_TYPE_VIRTUAL_CORES && !vcChanged){
            // We do not add statistics about workers reconfiguration
            // since they did not changed.
            ;
        }else{
            double ms = ticksToMilliseconds(getticks() - start,
                                            _p.archData.ticksPerNs);
            _reconfigurationStats.addSample(type, ms);
        }
    }
}

void Configuration::stopReconfigurationStatsTotal(ticks start){
    if(_p.statsReconfiguration){
        double ms = ticksToMilliseconds(getticks() - start,
                                        _p.archData.ticksPerNs);
        _reconfigurationStats.addSampleTotal(ms);
    }
}

void Configuration::setValues(const KnobsValues& values){
    bool vcChanged = virtualCoresWillChange(values);

    ticks totalStart = startReconfigurationStatsTotal();
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        ticks reconfigurationStart = startReconfigurationStatsKnob();

        // Start of the real reconfiguration
        if(values.areReal()){
            if(_knobs[i]->isLocked()){
                _knobs[i]->setRealValue(_knobs[i]->getRealValue());
            }else{
                _knobs[i]->setRealValue(values[(KnobType)i]);
            }
        }else if(values.areRelative()){
            if(_knobs[i]->isLocked()){
                _knobs[i]->setRelativeValue(0); // Only one value possible, any relative value would be ok.
            }else{
                _knobs[i]->setRelativeValue(values[(KnobType)i]);
            }
        }else{
            throw std::runtime_error("KnobsValues with undefined type.");
        }
        // End of the real reconfiguration
        stopReconfigurationStatsKnob(reconfigurationStart, (KnobType) i, vcChanged);
    }

    stopReconfigurationStatsTotal(totalStart);

    DEBUG("Changed knobs values.");
}

void Configuration::trigger(){
    for(size_t i = 0; i < TRIGGER_TYPE_NUM; i++){
        if(_triggers[i]){
            _triggers[i]->trigger();
        }
    }
}


ConfigurationExternal::ConfigurationExternal(const Parameters& p):
        Configuration(p){

    _knobs[KNOB_TYPE_VIRTUAL_CORES] = new KnobVirtualCores(p);
    _knobs[KNOB_TYPE_HYPERTHREADING] = new KnobHyperThreading(p);
    _knobs[KNOB_TYPE_MAPPING] = new KnobMappingExternal(p,
                                                *((KnobVirtualCores*) _knobs[KNOB_TYPE_VIRTUAL_CORES]),
                                                *((KnobHyperThreading*) _knobs[KNOB_TYPE_HYPERTHREADING]));
    _knobs[KNOB_TYPE_FREQUENCY] = new KnobFrequency(p,
                                                    *((KnobMappingExternal*) _knobs[KNOB_TYPE_MAPPING]));

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = NULL;
}


ConfigurationFarm::ConfigurationFarm(const Parameters& p,
                                     Smoother<MonitoredSample> const* samples,
                                     AdaptiveNode* emitter,
                                     std::vector<AdaptiveNode*> workers,
                                     AdaptiveNode* collector,
                                     ff::ff_gatherer* gt,
                                     volatile bool* terminated):
        Configuration(p), _numServiceNodes(0){

    _knobs[KNOB_TYPE_VIRTUAL_CORES] = new KnobVirtualCoresFarm(p,
                                          emitter, collector, gt, workers,
                                          terminated);

    _knobs[KNOB_TYPE_HYPERTHREADING] = new KnobHyperThreading(p);
    _knobs[KNOB_TYPE_MAPPING] = new KnobMappingFarm(p,
                                                *((KnobVirtualCoresFarm*) _knobs[KNOB_TYPE_VIRTUAL_CORES]),
                                                *((KnobHyperThreading*) _knobs[KNOB_TYPE_HYPERTHREADING]),
                                                emitter, collector);
    _knobs[KNOB_TYPE_FREQUENCY] = new KnobFrequency(p,
                                                    *((KnobMappingFarm*) _knobs[KNOB_TYPE_MAPPING]));

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = new TriggerQBlocking(p.triggerQBlocking,
                                                              p.thresholdQBlocking,
                                                              samples,
                                                              emitter);

    if(emitter){++_numServiceNodes;}
    if(collector){++_numServiceNodes;}
}


KnobsValues getRealValues(const Configuration& configuration, const KnobsValues& values){
    KnobsValues real(KNOB_VALUE_REAL);

    if(values.areRelative()){
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            double realv;
            assert(configuration.getKnob((KnobType)i)->getRealFromRelative(values[(KnobType)i], realv));
            real[(KnobType)i] = realv;
        }
    }else{
        real = values;
    }
    return real;
}

std::vector<AdaptiveNode*> convertWorkers(ff::svector<ff::ff_node*> w){
    std::vector<AdaptiveNode*> r;
    for(size_t i = 0; i < w.size(); i++){
        r.push_back(dynamic_cast<AdaptiveNode*>(w[i]));
    }
    return r;
}

}
