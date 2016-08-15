/*
 * manager-multi.cpp
 *
 * Created on: 22/07/2016
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

#include "manager-multi.hpp"
#include "external/Mammut/mammut/mammut.hpp"

#include <map>

using namespace std;
using namespace mammut;
using namespace mammut::topology;
using namespace mammut::cpufreq;

/**
 * It works only if we have htlevel == 1.
 */

#ifdef DEBUG_MANAGER_MULTI
#undef DEBUG
#define DEBUG(x) do { cerr << "[ManagerMulti] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#undef DEBUG
#define DEBUG(x)
#define DEBUGB(x)
#endif


namespace nornir{

template<typename A, typename B>
std::pair<B, A> flipPair(const std::pair<A,B> &p)
{
    return std::pair<B, A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B, A> flipMap(const std::map<A,B> &src)
{
    std::multimap<B, A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flipPair<A,B>);
    return dst;
}

ManagerMulti::ManagerMulti():_qIn(10), _qOut(10){
    _qIn.init();
    _qOut.init();
    _topology = _m.getInstanceTopology();
    _cpufreq = _m.getInstanceCpuFreq();
    vector<PhysicalCore*> allPc = _topology->getPhysicalCores();
    for(size_t i = 0; i < allPc.size(); i++){
        _allCores.push_back(allPc.at(i)->getPhysicalCoreId());
    }
}

void ManagerMulti::addManager(Manager* m){
    // Send the new manager
    while(!_qIn.push(m)){;}
}

Manager* ManagerMulti::getTerminatedManager(){
    Manager *m;
    if(_qOut.pop((void**) &m)){
        return m;
    }else{
        return NULL;
    }
}

void ManagerMulti::inhibitAll(Manager* except){
    for(size_t i = 0; i < _activeManagers.size(); i++){
        Manager* currentManager = _activeManagers.at(i);
        if(currentManager != except){
            currentManager->inhibit();
        }
    }
}

void ManagerMulti::disinhibitAll(){
    for(size_t i = 0; i < _activeManagers.size(); i++){
        _activeManagers.at(i)->disinhibit();
    }
}

vector<PhysicalCoreId> ManagerMulti::getAvailablePhysicalCores() const{
    vector<PhysicalCoreId> availablePc;
    for(size_t i = 0; i < _allCores.size(); i++){
        bool alreadyAllocated = false;
        for(auto it = _allocatedCores.begin(); it != _allocatedCores.end(); it++){
            if(utils::contains(it->second, _allCores.at(i))){
                alreadyAllocated = true;
                break;
            }
        }
        if(!alreadyAllocated){
            availablePc.push_back(_allCores.at(i));
        }
    }
    return availablePc;
}

map<KnobsValues, double> ManagerMulti::invertMap(const map<KnobsValues, double>& map) const{
    std::map<KnobsValues, double> mapr;
    for(auto it = map.begin(); it != map.end(); it++){
        mapr[it->first] = 1.0 / it->second;
    }
    return mapr;
}

void ManagerMulti::updateAllocations(Manager* m){
    std::map<KnobsValues, double> primaryValues;
    std::map<KnobsValues, double> secondaryValues;

    /**
     *  We want lower=better for both metrics. 
     *    For power consumption it is already true. 
     *    For bandwidth we need to invert it since it is higher=better.
     **/
    double primaryBound = 0.0;
    switch(m->_p.contractType){
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_BANDWIDTH:{
            primaryBound = 1.0 / m->_p.requiredBandwidth;
            primaryValues = invertMap(((SelectorPredictive*) m->_selector)->getPrimaryPredictions());
            secondaryValues = ((SelectorPredictive*) m->_selector)->getSecondaryPredictions();
        }break;
        case CONTRACT_POWER_BUDGET:{
            primaryBound = m->_p.powerBudget;
            secondaryValues = invertMap(((SelectorPredictive*) m->_selector)->getSecondaryPredictions());
            primaryValues = ((SelectorPredictive*) m->_selector)->getPrimaryPredictions();
        }break;
        default:{
            throw std::runtime_error("updateAllocations: Contract type not supported.");
        }break;
    }

    std::multimap<double, KnobsValues> unfeasible;
    std::multimap<double, KnobsValues> sortedSecondary;
    for(auto it = primaryValues.begin(); it != primaryValues.end(); ){
        KnobsValues kv = it->first;
        if(it->second <= primaryBound){
            // Insert the corresponding entry in secondaryValues
            sortedSecondary.insert(std::pair<double, KnobsValues>(secondaryValues.at(kv), kv));
            // Delete the corresponding entry in primaryValues
            auto del = it;
            it++;
            primaryValues.erase(del);
        }else{
            // Insert unfeasible solutions according to their distance from the bound
            unfeasible.insert(std::pair<double, KnobsValues>(it->second - primaryBound, kv));
            it++;
        }
    }
    std::vector<KnobsValues> allocation;
    for(auto it = sortedSecondary.begin(); it != sortedSecondary.end();  it++){
        allocation.push_back(it->second);
    }
    for(auto it = unfeasible.begin(); it != unfeasible.end();  it++){
        allocation.push_back(it->second);
    }
    _allocations[m] = allocation;
}

void ManagerMulti::updateAllocations(){
    for(size_t i = 0; i < _activeManagers.size(); i++){
        updateAllocations(_activeManagers[i]);
    }
    std::vector<std::vector<size_t>> values;
    std::vector<size_t> accum;
    for(size_t i = 0; i < _allocations.size(); i++){
        std::vector<size_t> tmp;
        for(size_t j = 0; j < _allocations.begin()->second.size(); j++){
            tmp.push_back(j);
        }
        values.push_back(tmp);
    }
    _allocationsCombinations.clear();
    combinations(values, 0, accum);
}

void ManagerMulti::combinations(std::vector<std::vector<size_t> > array, size_t i, std::vector<size_t> accum){
    if(i == array.size()){
        _allocationsCombinations.push_back(accum);
    }else{
        std::vector<size_t> row = array.at(i);
        for(size_t j = 0; j < row.size(); ++j){
            std::vector<size_t> tmp(accum);
            tmp.push_back(row[j]);
            combinations(array, i+1, tmp);
        }
    }
}

std::vector<size_t> ManagerMulti::findBestAllocation(){
    updateAllocations();
    std::vector<size_t> best;
    size_t minWeight = std::numeric_limits<size_t>::max();
    for(size_t i = 0; i < _allocationsCombinations.size(); i++){
        std::vector<size_t> allocation = _allocationsCombinations.at(i);
        size_t pos = 0;
        if(!allocation.size()){return best;}
        Frequency previousFreq = 0;
        size_t numPhysicalCores = 0;
        size_t weight = 0;
        bool validAllocation = true;
        for(auto it = _allocations.begin(); it != _allocations.end(); it++){
            weight += allocation.at(pos);
            KnobsValues kv = it->second.at(allocation.at(pos));
            size_t numCores;
            Frequency currentFreq;
            if(kv.areRelative()){
                double tmp;
                assert(it->first->_configuration->getKnob(KNOB_TYPE_FREQUENCY)->getRealFromRelative(kv[KNOB_TYPE_FREQUENCY], tmp));
                currentFreq = tmp;
                assert(it->first->_configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealFromRelative(kv[KNOB_TYPE_VIRTUAL_CORES], tmp));
                numCores = tmp;
            }else{
                currentFreq = kv[KNOB_TYPE_FREQUENCY];
                numCores = kv[KNOB_TYPE_VIRTUAL_CORES];
            }
            // Only keep combinations on the same frequency
            if(previousFreq && currentFreq != previousFreq){
                validAllocation = false;
                break;
            }
            previousFreq = currentFreq;
            numPhysicalCores += numCores;
            numPhysicalCores += it->first->_configuration->getNumServiceNodes();
            ++pos;
        }
        if(numPhysicalCores > _allCores.size()){
            validAllocation = false;
        }
        // We want the configuration such that the sum of the indexes
        // is minimal.
        if(validAllocation && weight < minWeight){
            minWeight = weight;
            best = allocation;
        }
    }
    if(minWeight == std::numeric_limits<size_t>::max()){
        throw std::runtime_error("No valid allocations found.");
    }
    return best;
}

void ManagerMulti::applyNewAllocation(){
    std::vector<size_t> alloc = findBestAllocation();
    DEBUG("Best allocation found: " << alloc);
    size_t pos = 0;
    size_t nextCoreId = 0;
    Manager* man;
    for(auto it = _allocations.begin(); it != _allocations.end(); it++){
        man = it->first;
        KnobsValues kv = it->second.at(alloc.at(pos));
        DEBUG("Allocations: " << it->second);
        DEBUG("Allocation: " << man << " " << kv << " " << ((SelectorPredictive*)man->_selector)->getPrimaryPrediction(kv) << " " << ((SelectorPredictive*)man->_selector)->getSecondaryPrediction(kv));
        size_t numCores;
        if(kv.areRelative()){
            double tmp;
            assert(man->_configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealFromRelative(kv[KNOB_TYPE_VIRTUAL_CORES], tmp));
            numCores = tmp;
        }else{
            numCores = kv[KNOB_TYPE_VIRTUAL_CORES];
        }
        numCores += man->_configuration->getNumServiceNodes();

        vector<PhysicalCoreId> cores;
        cores.reserve(numCores);
        for(size_t i = 0; i < numCores; i++){
            cores.push_back(_allCores.at(i + nextCoreId));
        }
        DEBUG("Forcing manager " << man << " to " << cores);
        man->allowPhysicalCores(cores);
        //man->_selector->forceConfiguration(kv);
        man->act(kv, true);
        ((SelectorPredictive*) man->_selector)->updatePredictions(kv);
        ++pos;
        nextCoreId += numCores;
    }
}

void ManagerMulti::applyWattsCorrection(){
    for(size_t i = 0; i < _activeManagers.size(); i++){
        _activeManagers.at(i)->updateWattsCorrection();
    }
}

PredictorLinearRegression* ManagerMulti::getPowerPredictor(Manager* m) const{
    switch(m->_p.contractType){
        case CONTRACT_PERF_BANDWIDTH:
        case CONTRACT_PERF_COMPLETION_TIME:
        case CONTRACT_PERF_UTILIZATION:{
            return ((PredictorLinearRegression*) ((SelectorPredictive*) m->_selector)->getSecondaryPredictor());
        }break;
        case CONTRACT_POWER_BUDGET:{
            return ((PredictorLinearRegression*) ((SelectorPredictive*) m->_selector)->getPrimaryPredictor());
        }break;
        default:{
            throw std::runtime_error("Unknown contract type.");
        }
    }
}

double ManagerMulti::getInactivePowerParameter(Manager* m) const{
    return getPowerPredictor(m)->getInactivePowerParameter();
}

void ManagerMulti::run(){
    while(true){
        Manager* m;
        if(_qIn.pop((void**) &m)){
            DEBUG("New manager arrived.");
            // New manager.
            /*
            assert(m->_p.contractType == CONTRACT_PERF_BANDWIDTH ||
                   m->_p.contractType == CONTRACT_PERF_COMPLETION_TIME ||
                   m->_p.contractType == CONTRACT_PERF_UTILIZATION);
            */
            assert(m->_p.strategySelection == STRATEGY_SELECTION_LEARNING);
            assert(m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_AMDAHL ||
                    m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_AMDAHL_MAPPING ||
                    m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL ||
                    m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL_MAPPING ||
                    m->_p.strategyPredictionPower == STRATEGY_PREDICTION_POWER_LINEAR ||
                    m->_p.strategyPredictionPower == STRATEGY_PREDICTION_POWER_LINEAR_MAPPING);

            _activeManagers.push_back(m);
            _allocatedCores[m] = vector<PhysicalCoreId>();
            if(_activeManagers.size() > 1){
                inhibitAll(m); //TODO: Inhibition should freeze the application.
            }
            DEBUG("All the other managers have been inhibited.");
            vector<PhysicalCoreId> cores = getAvailablePhysicalCores();
            DEBUG("New manager (" << m << ") can run on: " << cores);
            m->allowPhysicalCores(cores);
            m->_selector->setCalibrationCoordination();
            m->_selector->allowCalibration();
            m->start();
            DEBUG("Manager started.");
            // Wait for calibration termination.
            while(m->_selector->isCalibrating() || !m->_selector->getTotalCalibrationTime()){;}
            DEBUG("Calibration terminated.");

            DEBUG(m << " Inactive power parameter: " << getInactivePowerParameter(m));

            if(_activeManagers.size() > 1){
                applyNewAllocation();
                DEBUG("Best allocation applied.");
                // TODO Add external contributions to perf models.
                //applyWattsCorrection();
                // Disinhibit all the managers
                disinhibitAll();
                DEBUG("All managers disinhibited.");
            }
        }else{
            _allocatedCores.clear();
            // Manage already present managers.
            for(size_t i = 0; i < _activeManagers.size(); ){
                Manager* m = _activeManagers.at(i);
                if(!m->running()){
                    // Manager terminated
                    std::vector<PhysicalCoreId> released = m->getUsedCores();
                    _allocatedCores.erase(m);
                    _activeManagers.erase(_activeManagers.begin() + i);
                    if(_activeManagers.size()){
                        //applyWattsCorrection(); //TODO Possible contention on samples
                    }
                    while(!_qOut.push((void*) m)){;}
                }else{
                    i++;
                    if(m->_selector->isCalibrating()){
                        DEBUG("Manager " << m << " wants to calibrate.");
                        if(_activeManagers.size() > 1){
                            inhibitAll(m);
                        }
                        m->_selector->allowCalibration();
                        while(m->_selector->isCalibrating()){;}
                        DEBUG("Calibration terminated.");
                        DEBUG(m << " Inactive power parameter: " << getInactivePowerParameter(m));
                        if(_activeManagers.size() > 1){
                            applyNewAllocation();
                            DEBUG("Best allocation applied.");
                            // TODO Add external contributions to perf models.
                            //applyWattsCorrection();
                            // Disinhibit all the other managers
                            disinhibitAll();
                            DEBUG("All managers disinhibited.");
                        }
                    }else{
                        std::vector<PhysicalCoreId> pc;
                        do{
                            pc = m->getUsedCores();
                        }while(pc.size() == 0);
                        DEBUG("Manager " << m << " uses cores: " << pc);
                        _allocatedCores[m] = pc;
                        m->allowPhysicalCores(pc);
                    }
                }
            }
        }
        sleep(1);
    }
}

} // End namespace
