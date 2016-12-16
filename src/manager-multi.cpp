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
#include <signal.h>
#include <sys/types.h>

using namespace std;
using namespace mammut;
using namespace mammut::topology;
using namespace mammut::cpufreq;

/**
 * It works only if we have htlevel == 1.
 */

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_MANAGER_MULTI
#define DEBUG(x) do { cerr << "[ManagerMulti] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
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

ManagerMulti::ManagerMulti(ManagerMultiConfiguration configuration):
        _configuration(configuration), _qIn(10), _qOut(10),
        _power(new MovingAverageSimple<double>(MAX_POWER_VIOLATION_SECONDS)){
    if(_configuration.powerCap < 0){
        throw std::runtime_error("[ManagerMulti]: powerCap must be >= 0.");
    }
    _qIn.init();
    _qOut.init();
    _topology = _m.getInstanceTopology();

    if(configuration.useVirtualCores){
        vector<VirtualCore*> allCores = _topology->getVirtualCores();
        for(size_t i = 0; i < allCores.size(); i++){
            _allCores.push_back(allCores.at(i)->getVirtualCoreId());
        }
    }else{
        vector<PhysicalCore*> allCores = _topology->getPhysicalCores();
        for(size_t i = 0; i < allCores.size(); i++){
            _allCores.push_back(allCores.at(i)->getVirtualCore()->getVirtualCoreId());
        }
    }
}

ManagerMulti::~ManagerMulti(){
    delete _power;
}

void ManagerMulti::addManager(Manager* m, double minPerformanceRequired){
    // Send the new manager
    checkManagerSupported(m);
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

void ManagerMulti::allowCalibration(Manager* m){
    DEBUG("Manager (" << m << ") wants to calibrate.");
    vector<VirtualCoreId> cores = getAvailableCores();
    m->allowCores(cores);
    DEBUG("Manager (" << m << ") can calibrate on: " << cores);
    m->_selector->allowCalibration();
}

void ManagerMulti::waitForCalibration(Manager* m){
    // We need to check the total calibration time
    // to be sure that the calibration at least
    // started the first time.
    while(m->_selector->isCalibrating() ||
          !m->_selector->getTotalCalibrationTime()){
        ;
    }
    DEBUG("Manager (" << m << ") terminated its calibration.");
}

void ManagerMulti::inhibitAll(Manager* except){
    for(auto it : _managerData){
        Manager* currentManager = it.first;
        if(currentManager != except){
            currentManager->inhibit();
        }
    }
    DEBUG("Inhibition done.");
}

void ManagerMulti::disinhibitAll(){
    for(auto it : _managerData){
        Manager* currentManager = it.first;
        currentManager->disinhibit();
    }
    DEBUG("All managers disinhibited.");
}

void ManagerMulti::shrinkAll(Manager* except){
    for(auto it : _managerData){
        Manager* currentManager = it.first;
        if(currentManager != except){
            currentManager->shrink(_configuration.shrink);
            switch(_configuration.shrink){
                case CALIBRATION_SHRINK_AGGREGATE:{
                    it.second.allocatedCores.clear();
                    // TODO: At the moment I put everything on the same core.
                    // It would maybe be better
                    // to have each different application on a different core.
                    it.second.allocatedCores.push_back(_topology->getVirtualCores().back()->getVirtualCoreId());
                }break;
                case CALIBRATION_SHRINK_PAUSE:{
                    it.second.allocatedCores.clear();
                }break;
                default:{
                    ;
                }break;
            }
        }
    }
    DEBUG("Shrinking done.");
}

vector<VirtualCoreId> ManagerMulti::getAvailableCores() const{
    vector<VirtualCoreId> availableCores;
    for(size_t i = 0; i < _allCores.size(); i++){
        bool alreadyAllocated = false;
        for(auto it : _managerData){
            if(utils::contains(it.second.allocatedCores, _allCores.at(i))){
                alreadyAllocated = true;
                break;
            }
        }
        if(!alreadyAllocated){
            availableCores.push_back(_allCores.at(i));
        }
    }
    return availableCores;
}

void ManagerMulti::updateAllocations(Manager* m){
    std::map<KnobsValues, double> primaryValues;
    std::map<KnobsValues, double> secondaryValues;

    double primaryBound = 0.0;
    primaryBound = m->_p.requiredBandwidth;
    primaryValues = ((SelectorPredictive*) m->_selector)->getPrimaryPredictions();
    secondaryValues = ((SelectorPredictive*) m->_selector)->getSecondaryPredictions();

    ManagerData md = _managerData[m];
    double referencePerformance;
    if(m->_p.contractType == CONTRACT_PERF_MAX){
        double maxPerformance = 0.0;
        for(auto it : primaryValues){
            if(it.second > maxPerformance){
                maxPerformance = it.second;
            }
        }
        referencePerformance = maxPerformance;
    }else{
        referencePerformance = m->_p.requiredBandwidth;
    }
    md.minPerf = (md.minPerfReqPerc/100.0) * referencePerformance;

    std::multimap<double, KnobsValues> unfeasible;
    std::multimap<double, KnobsValues> sortedSecondary;
    for(auto it = primaryValues.begin(); it != primaryValues.end(); ){
        KnobsValues kv = it->first;
        if(it->second >= primaryBound){
            // Insert the corresponding entry in secondaryValues
            // since is a map, they will be kept sorted from the lower power consuming
            // to the higher power consuming.
            sortedSecondary.insert(std::pair<double, KnobsValues>
                                  (secondaryValues.at(kv), kv));
            // Delete the corresponding entry in primaryValues
            it = primaryValues.erase(it);
        }else{
            // Insert unfeasible solutions according to their relative performance in
            // percentage (from lowest to highest).
            double relativePerf = (it->second / referencePerformance) * 100;
            unfeasible.insert(std::pair<double, KnobsValues>(relativePerf, kv));
            it++;
        }
    }
    std::vector<std::pair<KnobsValues, double> > allocations;
    // First insert the solutions that satisfies the primary bound (sorted from 
    // the best to the worst secondary value).
    for(auto it = sortedSecondary.begin(); it != sortedSecondary.end();  it++){
        // Here all these solutions satisfy the constraint, and we only want to
        // order according to the secondary value. So we are not interested on
        // the relative performance. So we put 100%.
        double relativePerf = 100;
        allocations.push_back(std::pair<KnobsValues, double>(it->second, relativePerf));
    }
    // Then we insert the unfeasible solutions (i.e. that violate the primary bound)
    // sorted from the most performing to the least performing.
    // We scan on the reverse direction since they are ordered from the least
    // to the most performing.
    for(auto it = unfeasible.rbegin(); it != unfeasible.rend(); it++){
        double relativePerf = (it->first / referencePerformance) * 100;
        allocations.push_back(std::pair<KnobsValues, double>(it->second, relativePerf));
    }
    _managerData[m].allocations = allocations;
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

double ManagerMulti::getQuality(const std::vector<size_t>& indexes)const {
    switch(_configuration.qualityEstimation){
        case QUALITY_ESTIMATION_PERFORMANCE:{
            // For each manager, we check how much performance degradation (in percentage)
            // we would get. Then, we return the average between these degradations.
            double avgRelPerf = 0;
            size_t pos = 0;
            for(auto it : _managerData){
                avgRelPerf += it.second.allocations.at(indexes.at(pos)).second;
                ++pos;
            }
            avgRelPerf /= _managerData.size();
            return avgRelPerf;
        }break;
        case QUALITY_ESTIMATION_PREFERENCE:{
            // At the moment we are considering as quality metric the sum of
            // the positions of the chosen configurations in the corresponding priority lists
            // (the lower the better).
            // If we have 3 managers and we select the second choice for manager A, the
            // tenth choice for manager B and the fifth choice for manager C, then this
            // specific allocation will have a quality of 2+10+5 = 17. Accordingly, in
            // this specific case the best allocation would be the one with weight 3, i.e.
            // the one in which to each manager we assign its first (most preferred) choice.
            // We could also consider other metrics (e.g. average).
            double r = 0;
            for(size_t x : indexes){r += x;}
            // Since it is the lowest the better but we want the higher the better.
            return std::numeric_limits<double>::max() - r;
        }break;
        default:{
            throw std::runtime_error("Unsupported quality estimation algorithm.");
        }
    }
}

double ManagerMulti::estimatePower(const std::vector<size_t>& indexes) const{
    return 0; //TODO
}

void ManagerMulti::calibrate(Manager* m, bool start){
    inhibitAll(m);
    shrinkAll(m);
    allowCalibration(m);
    if(start){
        m->start(); DEBUG("Manager started.");
    }
    waitForCalibration(m);
    applyNewAllocation();
    m->inhibit(); DEBUG(m << " inhibited.");
    for(auto it : _managerData){
        if(it.first != m){
            it.first->stretch(_configuration.shrink);
        }
    }
    DEBUG("Everyone stretched.");
    sleep(1);
    DEBUG("Ready to update the models.");
    for(auto it : _managerData){
        it.first->updateModelsInterference();
        it.first->waitModelsInterferenceUpdate();
    }
    DEBUG("Models updated.");
    disinhibitAll();
}

std::vector<size_t> ManagerMulti::findBestAllocation(){
    updateAllocations();
    Manager* currentManager;
    KnobsValues kv;
    Frequency currentFreq;
    // Feasible solutions, sorted from the best to the worst
    std::multimap<double, std::vector<size_t> > feasibleSolutions;
    for(std::vector<size_t> indexes : _allocationsCombinations){
        size_t pos = 0;
        if(!indexes.size()){
            throw std::runtime_error("FATAL ERROR: No indexes.");
        }
        Frequency previousFreq = 0;
        size_t numCores = 0;
        bool validAllocation = true;
        std::pair<KnobsValues, double> allocation;
        size_t allocationPosition;
        for(auto it : _managerData){
            currentManager = it.first;
            const ManagerData& currentManagerData = it.second;
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
            numCores += numCores;
            numCores += currentManager->_configuration->getNumServiceNodes();
            ++pos;
        }
        if(numCores > _allCores.size() ||
           estimatePower(indexes) > _configuration.powerCap){
            validAllocation = false;
        }

        if(validAllocation){
            DEBUG("Allocation " << indexes << " has quality " << getQuality(indexes));
            feasibleSolutions.insert(std::pair<double, std::vector<size_t> >
                (getQuality(indexes), indexes));
        }else{
            DEBUG("Allocation " << indexes << " is not valid.");
            // Just to avoid having an empty map where there are
            // no feasible solutions. If it is not valid, we consider
            // the quality as negative.
            feasibleSolutions.insert(std::pair<double, std::vector<size_t> >
                (-1, indexes));
        }
    }

    // First try to find a solution that doesn't violate
    // any additional requirement. If does not exists,
    // just return the one with maximum quality.
    // We need to scan in the reverse direction since they are ordered from
    // the one with lowest quality to the one with highest quality.
    std::pair<KnobsValues, double> allocation;
    for(auto indexes = feasibleSolutions.rbegin();
             indexes != feasibleSolutions.rend();
             indexes++){
        size_t pos = 0;
        size_t allocationPosition;
        bool feasible = true;
        for(auto it : _managerData){
            currentManager = it.first;
            const ManagerData& currentManagerData = it.second;
            allocationPosition = indexes->second.at(pos);
            allocation = currentManagerData.allocations.at(allocationPosition);
            // Check that we still satisfy the additional requirement
            // set on the global manager.
            if(allocation.second < currentManagerData.minPerf){
                feasible = false;
                break;
            }
            ++pos;
        }
        // Since we are iterating from the best to the worst,
        // as soon as we find a feasible one we return.
        if(feasible){return indexes->second;}
    }

    DEBUG("No feasible solutions found.");
    // If we are here, there are no solutions that satisfy the
    // additional performance constraints, so we return the
    // one with maximum quality.
    return feasibleSolutions.end()->second;
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
        //DEBUG("Allocations: " << it.second.allocations);
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

        vector<VirtualCoreId> cores;
        cores.reserve(numCores);
        for(size_t i = 0; i < numCores; i++){
            cores.push_back(_allCores.at(i + nextCoreId));
        }
        DEBUG("Forcing manager " << man << " to " << cores);
        man->allowCores(cores);
        //man->_selector->forceConfiguration(kv);
        man->act(kv, true);
        ((SelectorPredictive*) man->_selector)->updatePredictions(kv);
        ++pos;
        nextCoreId += numCores;
    }
    DEBUG("Best allocation applied.");
}

void ManagerMulti::checkManagerSupported(Manager* m){
    // It is meaningless to control power for individual applications if
    // other applications are running on the system and want to be
    // controlled as well.
    assert(m->_p.contractType != CONTRACT_POWER_BUDGET);

    assert(m->_p.strategySelection == STRATEGY_SELECTION_LEARNING);
    assert(m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_AMDAHL ||
            m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL ||
            m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP);
}

void ManagerMulti::run(){
    mammut::energy::Counter* joulesCounter = _m.getInstanceEnergy()->getCounter();
    double lastSampleTime = mammut::utils::getMillisecondsTime();
    double lastJoules = joulesCounter->getJoules();
    double currentSampleTime, currentJoules, currentWatts;
    while(true){
        SubmittedManager* sm;
        if(_qIn.pop((void**) &sm)){
            ManagerData md;
            Manager* m = sm->manager;
            DEBUG("Manager (" << m << ") arrived with contract " << m->_p.contractType);
            md.minPerfReqPerc = sm->minPerf;
            _managerData[m] = md;
            m->_selector->setCalibrationCoordination();
            calibrate(m, true);
            delete sm;
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
                        calibrate(m);
                    }else{
                        std::vector<VirtualCoreId> pc;
                        do{
                            pc = m->getUsedCores();
                        }while(pc.size() == 0);
                        DEBUG("Manager " << m << " uses cores: " << pc);
                        _managerData[m].allocatedCores = pc;
                        m->allowCores(pc);
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
        if(_power->average() > _configuration.powerCap){
            DEBUG("Cap violated (" << _power->average() << ">" << _configuration.powerCap << ". "
                  "Falling back to RAP.");
            ; //TODO Fallback to RAPL
        }
        sleep(1);
    }
}

} // End namespace
