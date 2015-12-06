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

namespace adpff{

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

}
