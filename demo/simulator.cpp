/*
 * simulator.hpp
 *
 * Created on: 26/05/2016
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


/**
 * Basic test for the simulator.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "../src/interface.hpp"

using namespace nornir;

volatile bool terminationFlag;

class DummyEmitter: public AdaptiveNode {
public:
    DummyEmitter(){;}

    void* svc(void *) {
        if(terminationFlag){
            return NULL;
        }else{
            return (void*) 0x1;
        }
    }
};

class DummyNode: public AdaptiveNode {
public:
    DummyNode(){;}

    void* svc(void *) {return (void*) 0x1;}
};

void updateParameters(Parameters& p, double requirementA, double requirementB = 0){
    switch(p.contractType){
        case CONTRACT_PERF_UTILIZATION:{
            p.underloadThresholdFarm = requirementA;
            p.overloadThresholdFarm = requirementB;
        }break;
        case CONTRACT_PERF_BANDWIDTH:{
            p.requiredBandwidth = requirementA;
        }break;
        case CONTRACT_PERF_COMPLETION_TIME:{
            p.requiredCompletionTime = requirementA;
        }break;
        case CONTRACT_POWER_BUDGET:{
            p.powerBudget = requirementA;
        }break;
        default:{
            ;
        }break;
    }
}

std::vector<double> getPowerBounds(std::string confData, mammut::Mammut& mammut){
    std::vector<double> r;
    std::vector<std::string> lines = mammut::utils::readFile(confData);
    double minPw = std::numeric_limits<double>::max(),
           maxPw = 0, pw;
    for(size_t i = 1; i < lines.size(); i++){
        // Check if the frequency in the file is a supported frequency. (e.g. tu avoid turbo boost frequencies.)
        if(!mammut::utils::contains(mammut.getInstanceCpuFreq()->getDomains().back()->getAvailableFrequencies(),
                                    (Frequency) atof(mammut::utils::split(lines.at(i), '\t')[1].c_str()))){
            continue;
        }
        pw = atof(mammut::utils::split(lines.at(i), '\t')[4].c_str());
        if(pw < minPw){
            minPw = pw;
        }
        if(pw > maxPw){
            maxPw = pw;
        }
    }

    for(size_t i = 10; i <= 100; i += 20){
        r.push_back(minPw + (maxPw - minPw)*(i/100.0));
    }
    return r;
}

std::vector<double> getBandwidthBounds(std::string confData, mammut::Mammut& mammut){
    std::vector<double> r;
    std::vector<std::string> lines = mammut::utils::readFile(confData);
    double minBw = std::numeric_limits<double>::max(),
           maxBw = 0, bw;
    for(size_t i = 1; i < lines.size(); i++){
        // Check if the frequency in the file is a supported frequency. (e.g. tu avoid turbo boost frequencies.)
        if(!mammut::utils::contains(mammut.getInstanceCpuFreq()->getDomains().back()->getAvailableFrequencies(),
                                    (Frequency) atof(mammut::utils::split(lines.at(i), '\t')[1].c_str()))){
            continue;
        }
        bw = 1.0 / atof(mammut::utils::split(lines.at(i), '\t')[2].c_str());
        if(bw < minBw){
            minBw = bw;
        }
        if(bw > maxBw){
            maxBw = bw;
        }
    }

    for(size_t i = 10; i <= 100; i += 20){
        r.push_back(minBw + (maxBw - minBw)*(i/100.0));
    }
    return r;
}

double getBest(std::string confData, ContractType contract, double bound){
    std::vector<std::string> lines = mammut::utils::readFile(confData);
    double bestBw = 0, bestPw = std::numeric_limits<double>::max(), bw, pw;
    double best = 0;
    for(size_t i = 1; i < lines.size(); i++){
        pw = atof(mammut::utils::split(lines.at(i), '\t')[4].c_str());
        bw = 1.0 / atof(mammut::utils::split(lines.at(i), '\t')[2].c_str());
        switch(contract){
            case CONTRACT_PERF_BANDWIDTH:{
                if(bw >= bound && pw < bestPw){
                    bestPw = pw;
                    best = bestPw;
                }
            }break;
            case CONTRACT_POWER_BUDGET:{
                if(pw <= bound && bw > bestBw){
                    bestBw = bw;
                    best = bestBw;
                }
            }break;
            default:{
                throw std::runtime_error("Unsupported contract.");
            }break;
        }
    }
    return best;
}

void buildFarm(ff_farm<>& farm, size_t numWorkers){
    farm.add_emitter(new DummyEmitter);
    std::vector<ff::ff_node*> w;
    for(size_t i = 0; i < numWorkers; ++i) w.push_back(new DummyNode);
    farm.add_workers(w);
    farm.add_collector(new DummyNode);
}

int main(int argc, char * argv[]) {
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " configurationFile [maxWorkers]" << std::endl;
        return -1;
    }

    std::string confData(argv[1]);

    uint numWorkers = 0;
    if(argc > 2){
        numWorkers = atoi(argv[2]);
    }else{
        std::vector<std::string> lines = mammut::utils::readFile(confData);
        for(size_t i = 1; i < lines.size(); i++){
            uint currNumWorkers = atoi(mammut::utils::split(lines.at(i), '\t')[0].c_str());
            if(currNumWorkers > numWorkers){
                numWorkers = currNumWorkers;
            }
        }
    }

    Parameters p("parameters.xml");
    std::vector<double> bounds;

    switch(p.contractType){
        case CONTRACT_PERF_BANDWIDTH:{
            bounds = getBandwidthBounds(confData, p.mammut);
        }break;
        case CONTRACT_POWER_BUDGET:{
            bounds = getPowerBounds(confData, p.mammut);
        }break;
        default:{
            throw std::runtime_error("Unsupported contract.");
        }break;
    }

    std::cout << "#Steps\tLoss" << std::endl;

    for(size_t i = 0; i < bounds.size(); i++){
        // Farm is not used, we just create it to let things work.
        ff_farm<> farm;
        buildFarm(farm, numWorkers);
        updateParameters(p, bounds.at(i));

        terminationFlag = false;
        ManagerFarm<> amf(&farm, p);
        SimulationResult sr = amf.simulate(std::string(argv[1]), &terminationFlag);
        //ATTENTION: No need to join. simulate() call will not start a new thread.
        double opt = getBest(confData, p.contractType, bounds.at(i));
        double loss = 0;
        bool sat = false;
        std::cout << sr.numSteps << "\t";

        switch(p.contractType){
            case CONTRACT_PERF_BANDWIDTH:{
                if(sr.currentBandwidth >= bounds.at(i)){
                    sat = true;
                    loss = (sr.currentPower - opt) / opt;
                }
            }break;
            case CONTRACT_POWER_BUDGET:{
                if(sr.currentPower <= bounds.at(i)){
                    sat = true;
                    loss = (opt - sr.currentBandwidth) / opt;
                }
            }break;
            default:{
                throw std::runtime_error("Unsupported contract.");
            }break;
        }

        if(sat){
            loss *= 100.0;
            std::cout << loss << std::endl;
        }else{
            std::cout << "MISS" << std::endl;
        }
    }

    return 0;
}

