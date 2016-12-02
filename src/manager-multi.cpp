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

class SubmittedManager{
public:
    Manager* manager;
    double minPerf;
    SubmittedManager(Manager* m, double minP):manager(m), minPerf(minP){;}
};

#define MAX_POWER_VIOLATION_SECONDS 10

ManagerMulti::ManagerMulti(double powerCap):
        _powerCap(powerCap), _qIn(10), _qOut(10),
        _power(new MovingAverageSimple<double>(MAX_POWER_VIOLATION_SECONDS)){
    if(_powerCap < 0){
        throw std::runtime_error("[ManagerMulti]: powerCap must be >= 0.");
    }
    _qIn.init();
    _qOut.init();
    _topology = _m.getInstanceTopology();
    vector<PhysicalCore*> allPc = _topology->getPhysicalCores();
    for(size_t i = 0; i < allPc.size(); i++){
        _allCores.push_back(allPc.at(i)->getPhysicalCoreId());
    }
}

ManagerMulti::~ManagerMulti(){
    delete _power;
}

void ManagerMulti::addManager(Manager* m, double minPerformanceRequired){
    // Send the new manager
    SubmittedManager* sm = new SubmittedManager(m, minPerformanceRequired);
    while(!_qIn.push(sm)){;}
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
    for(auto it : _managerData){
        Manager* currentManager = it.first;
        if(currentManager != except){
            currentManager->inhibit();
        }
    }
}

void ManagerMulti::disinhibitAll(){
    for(auto it : _managerData){
        it.first->disinhibit();
    }
}

vector<PhysicalCoreId> ManagerMulti::getAvailablePhysicalCores() const{
    vector<PhysicalCoreId> availablePc;
    for(size_t i = 0; i < _allCores.size(); i++){
        bool alreadyAllocated = false;
        for(auto it : _managerData){
            if(utils::contains(it.second.allocatedCores, _allCores.at(i))){
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

void ManagerMulti::updateAllocations(Manager* m){
    std::map<KnobsValues, double> primaryValues;
    std::map<KnobsValues, double> secondaryValues;

    double primaryBound = 0.0;
    double maxPerformance = 0.0;
    primaryBound = m->_p.requiredBandwidth;
    primaryValues = ((SelectorPredictive*) m->_selector)->getPrimaryPredictions();
    secondaryValues = ((SelectorPredictive*) m->_selector)->getSecondaryPredictions();
    for(auto it : primaryValues){
        if(it.second > maxPerformance){
            maxPerformance = it.second;
        }
    }
    ManagerData md = _managerData[m];
    if(m->_p.contractType == CONTRACT_PERF_MAX){
        md.minPerf = (md.minPerfReqPerc/100.0) * maxPerformance;
    }else{
        md.minPerf = (md.minPerfReqPerc/100.0) * m->_p.requiredBandwidth;
    }

    std::multimap<double, KnobsValues> unfeasible;
    std::multimap<double, KnobsValues> sortedSecondary;
    for(auto it = primaryValues.begin(); it != primaryValues.end(); ){
        KnobsValues kv = it->first;
        if(it->second <= primaryBound){
            // Insert the corresponding entry in secondaryValues
            sortedSecondary.insert(std::pair<double, KnobsValues>(secondaryValues.at(kv), kv));
            // Delete the corresponding entry in primaryValues
            it = primaryValues.erase(it);
        }else{
            // Insert unfeasible solutions according to their distance from the bound
            unfeasible.insert(std::pair<double, KnobsValues>(primaryBound - it->second, kv));
            it++;
        }
    }
    std::vector<std::pair<KnobsValues, double> > allocation;
    // First insert the solutions that satisfies the primary bound (sorted from 
    // the best to the worst secondary value).
    for(auto it = sortedSecondary.begin(); it != sortedSecondary.end();  it++){
        double perf = primaryValues[it->second];
        allocation.push_back(std::pair<KnobsValues, double>(it->second, perf));
    }
    // Then we insert the unfeasible solutions (i.e. that violate the primary bound)
    // sorted from the closer to the farther from the bound.
    for(auto it = unfeasible.begin(); it != unfeasible.end(); it++){
        double perf = primaryValues[it->second];
        allocation.push_back(std::pair<KnobsValues, double>(it->second, perf));
    }
    _managerData[m].allocations = allocation;
}

void ManagerMulti::updateAllocations(){
    for(auto it : _managerData){
        updateAllocations(it.first);
    }
    std::vector<std::vector<size_t>> values;
    std::vector<size_t> accum;
    for(auto it : _managerData){
        std::vector<size_t> tmp;
        for(size_t j = 0; j < it.second.allocations.size(); j++){
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

static double getQuality(const std::vector<size_t>& indexes){
    // At the moment we are considering as quality metric the sum of
    // the positions of the chosen configurations in the corresponding priority lists
    // (the lower the better).
    // If we have 3 managers and we select the second choice for manager A, the
    // tenth choice for manager B and the fifth choice for manager C, then this
    // specific allocation will have a quality of 2+10+5 = 17. Accordingly, in
    // this specific case the best allocation would be the one with weight 3, i.e.
    // the one in which to each manager we assign its first (most preferred) choice.
    // We could also consider other metrics (e.g. average).
    double r;
    for(size_t x : indexes){
        r += x;
    }
    return r;
}

std::vector<size_t> ManagerMulti::findBestAllocation(){
    updateAllocations();
    Manager* currentManager;
    const ManagerData& currentManagerData;
    KnobsValues kv;
    size_t numCores;
    Frequency currentFreq;
    // Feasible solutions, sorted from the best to the worst
    std::multimap<double, std::vector<size_t> > feasibleSolutions;
    for(std::vector<size_t> indexes : _allocationsCombinations){
        size_t pos = 0;
        if(!indexes.size()){
            throw std::runtime_error("FATAL ERROR: No indexes.");
        }
        Frequency previousFreq = 0;
        size_t numPhysicalCores = 0;
        bool validAllocation = true;
        std::pair<KnobsValues, double> allocation;
        size_t allocationPosition;
        for(auto it : _managerData){
            currentManager = it.first;
            currentManagerData = it.second;
            allocationPosition = indexes.at(pos);
            allocation = currentManagerData.allocations.at(allocationPosition);
            kv = allocation.first;
            if(kv.areRelative()){
                double tmp;
                assert(currentManager->_configuration->getKnob(KNOB_TYPE_FREQUENCY)->getRealFromRelative(kv[KNOB_TYPE_FREQUENCY], tmp));
                currentFreq = tmp;
                assert(currentManager->_configuration->getKnob(KNOB_TYPE_VIRTUAL_CORES)->getRealFromRelative(kv[KNOB_TYPE_VIRTUAL_CORES], tmp));
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
            numPhysicalCores += currentManager->_configuration->getNumServiceNodes();
            ++pos;
        }
        if(numPhysicalCores > _allCores.size()){
            validAllocation = false;
        }
        double weight = getQuality(indexes);
        if(validAllocation){
            feasibleSolutions[weight] = indexes;
        }else{
            // Just to avoid having an empty map where there are
            // no feasible solutions.
            feasibleSolutions[std::numeric_limits<double>::max()] = indexes;
        }
    }

    // First try to find a solution that doesn't violate
    // any additional requirement. If does not exists,
    // just return the one with minimum weight.
    for(auto indexes : feasibleSolutions){
        size_t pos = 0;
        size_t allocationPosition;
        std::pair<KnobsValues, double> allocation;
        bool feasible = true;
        for(auto it : _managerData){
            currentManager = it.first;
            currentManagerData = it.second;
            allocationPosition = indexes.second.at(pos);
            allocation = currentManagerData.allocations.at(allocationPosition);
            // Check that we still satisfy the additional requirement
            // set on the global manager.
            if(allocation.second < currentManagerData.minPerf){
                feasible = false;
                break;
            }
            ++pos;
        }
        // Since they are ordered from the best to the worst,
        // as soon as we find a feasible one we return.
        if(feasible){
            return indexes;
        }
    }

    DEBUG("No feasible solutions found.");
    // If we are here, there are no solutions that satisfy the
    // additional performance constraints, so we return the
    // one with minimum weight.
    return feasibleSolutions.begin()->second;
}

void ManagerMulti::applyNewAllocation(){
    std::vector<size_t> alloc = findBestAllocation();
    DEBUG("Best allocation found: " << alloc);
    size_t pos = 0;
    size_t nextCoreId = 0;
    Manager* man;
    for(auto it : _managerData){
        man = it.first;
        KnobsValues kv = it.second.allocations.at(alloc.at(pos)).first;
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

void ManagerMulti::run(){
    mammut::energy::Counter* joulesCounter = _m.getInstanceEnergy()->getCounter();
    double lastSampleTime = mammut::utils::getMillisecondsTime();
    double lastJoules = joulesCounter->getJoules();
    double currentSampleTime, currentJoules, currentWatts;
    while(true){
        Manager* m;
        SubmittedManager* sm;
        if(_qIn.pop((void**) &sm)){
            m = sm->manager;
            ManagerData md;
            md.minPerfReqPerc = sm->minPerf;
            delete sm;
            DEBUG("New manager arrived.");
            // New manager.
            // It is meaningless to control power for individual applications if
            // other applications are running on the system and want to be
            // controlled as well.
            assert(m->_p.contractType != CONTRACT_POWER_BUDGET);

            assert(m->_p.strategySelection == STRATEGY_SELECTION_LEARNING);
            assert(m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_AMDAHL ||
                    m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL ||
                    m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP ||
                    m->_p.strategyPredictionPower == STRATEGY_PREDICTION_POWER_LINEAR);

            _managerData[m] = md;
            if(_managerData.size() > 1){
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

            if(_managerData.size() > 1){
                applyNewAllocation();
                DEBUG("Best allocation applied.");
                // TODO Add external contributions to perf models.
                // Disinhibit all the managers
                disinhibitAll();
                DEBUG("All managers disinhibited.");
            }
        }else{
            // Manage already present managers.
            for(auto it = _managerData.begin(); it != _managerData.end() ; ){
                Manager* m = it->first;
                if(!m->running()){
                    // Manager terminated
                    it = _managerData.erase(it);
                    // TODO Remove external contributions from perf models.
                    while(!_qOut.push((void*) m)){;}
                }else{
                    ++it;
                    if(m->_selector->isCalibrating()){
                        DEBUG("Manager " << m << " wants to calibrate.");
                        if(_managerData.size() > 1){
                            inhibitAll(m);
                        }
                        m->_selector->allowCalibration();
                        while(m->_selector->isCalibrating()){;}
                        DEBUG("Calibration terminated.");
                        if(_managerData.size() > 1){
                            applyNewAllocation();
                            DEBUG("Best allocation applied.");
                            // TODO Add external contributions to perf models.
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
                        _managerData[m].allocatedCores = pc;
                        m->allowPhysicalCores(pc);
                    }
                }
            }
        }
        currentJoules = joulesCounter->getJoules();
        currentSampleTime = mammut::utils::getMillisecondsTime();
        currentWatts = (currentJoules - lastJoules)/((currentSampleTime - lastSampleTime)/1000.0);
        lastSampleTime = currentSampleTime;
        lastJoules = currentJoules;
        _power->add(currentWatts);
        // TODO In realtà bisognerebbe controllare che non lo sforiamo
        // per un periodo consecutivo sostenuto
        if(_power->average() > _powerCap){
            DEBUG("Cap violated (" << _power->average() << ">" << _powerCap << ". "
                  "Falling back to RAP.");
            ; //TODO Fallback to RAPL
        }
        sleep(1);
    }
}

} // End namespace
