/*
 * configuration.cpp
 *
 * Created on: 05/12/2015
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

#include "configuration.hpp"

#ifdef DEBUG_CONFIGURATION
#define DEBUG(x) do { cerr << "[Configuration] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace adpff{

FarmConfiguration::FarmConfiguration(const Parameters& p, AdaptiveNode* emitter,
        AdaptiveNode* collector, ff::ff_gatherer* gt,
        std::vector<AdaptiveNode*> workers,
        Smoother<MonitoredSample> const* samples):
        _p(p),
        _knobsChangeNeeded(false){
    /************************************************************/
    /*                           KNOBS                          */
    /************************************************************/
    _knobs[KNOB_TYPE_WORKERS] = new KnobWorkers(p.knobWorkers, emitter,
                                                collector, gt, workers);
    _knobs[KNOB_TYPE_MAPPING] = new KnobMapping(p.knobMapping,
                                                p.knobMappingEmitter,
                                                p.knobMappingCollector,
                                                p.knobHyperthreading,
                                                p.mammut,
                                                emitter,
                                                collector,
                                                *((KnobWorkers*)_knobs[KNOB_TYPE_WORKERS]));
    _knobs[KNOB_TYPE_FREQUENCY] = new KnobFrequency(p.knobFrequencies,
                                                    p.mammut,
                                                    p.turboBoost,
                                                    p.strategyUnusedVirtualCores,
                                                    *((KnobMapping*)_knobs[KNOB_TYPE_MAPPING]));

    std::vector<std::vector<double>> values;
    std::vector<double> accum;
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        values.push_back(_knobs[i]->getAllowedValues());
        if(_knobs[i]->needsCalibration()){
            _knobsChangeNeeded = true;
        }
    }
    combinations(values, 0, accum);

    /************************************************************/
    /*                         TRIGGERS                         */
    /************************************************************/
    _triggers[TRIGGER_TYPE_Q_BLOCKING] = new TriggerQBlocking(p.triggerQBlocking,
                                                              p.thresholdQBlocking,
                                                              samples,
                                                              emitter);
}

FarmConfiguration::~FarmConfiguration(){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        delete _knobs[i];
    }
}

//TODO: Works even if a std::vector is empty? (i.e. a knob has no values)
void FarmConfiguration::combinations(std::vector<std::vector<double> > array, size_t i, std::vector<double> accum){
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

bool FarmConfiguration::equal(KnobsValues values) const{
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

bool FarmConfiguration::knobsChangeNeeded() const{
    return _knobsChangeNeeded;
}

const std::vector<KnobsValues>& FarmConfiguration::getAllRealCombinations() const{
    return _combinations;
}

void FarmConfiguration::setFastReconfiguration(){
    if(_p.fastReconfiguration){
        ((KnobFrequency*) _knobs[KNOB_TYPE_FREQUENCY])->setRelativeValue(100.0);
    }
}

const Knob* FarmConfiguration::getKnob(KnobType t) const{
    return _knobs[t];
}

void FarmConfiguration::maxAllKnobs(){
    if(_p.contractType == CONTRACT_NONE || !_knobsChangeNeeded){
        return;
    }
    DEBUG("Maxing all the knobs.");
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[(KnobType) i]->setToMax();
    }
}

double FarmConfiguration::getRealValue(KnobType t) const{
    return _knobs[t]->getRealValue();
}

KnobsValues FarmConfiguration::getRealValues() const{
    KnobsValues kv(KNOB_VALUE_REAL);
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        kv[(KnobType) i] = getRealValue((KnobType) i);
    }
    return kv;
}

void FarmConfiguration::setRelativeValues(const KnobsValues& values){
    // Fast reconfiguration is valid only for knobs changed before
    // the frequency knob.
    assert(values.areRelative());
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[i]->setRelativeValue(values[(KnobType)i]);
    }
    DEBUG("Changed relative knobs values.");
}

void FarmConfiguration::setRealValues(const KnobsValues& values){
    // Fast reconfiguration is valid only for knobs changed before
    // the frequency knob.
    assert(values.areReal());
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _knobs[i]->setRealValue(values[(KnobType)i]);
    }
    DEBUG("Changed real knobs values.");
}

void FarmConfiguration::setValues(const KnobsValues& values){
    if(values.areReal()){
        setRealValues(values);
    }else if(values.areRelative()){
        setRelativeValues(values);
    }else{
        throw std::runtime_error("KnobsValues with undefined type.");
    }
}

void FarmConfiguration::trigger(){
    for(size_t i = 0; i < TRIGGER_TYPE_NUM; i++){
        _triggers[i]->trigger();
    }
}

}
