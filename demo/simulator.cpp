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

class DummyNode: public nornir::AdaptiveNode {
public:
    DummyNode(){;}

    void* svc(void *) {return (void*) 0x1;}
};

int main(int argc, char * argv[]) {
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " configurationFile" << std::endl;
    }

    uint numWorkers = 0;
    std::vector<std::string> lines = mammut::utils::readFile(std::string(argv[1]));
    for(size_t i = 0; i < lines.size(); i++){
        uint currNumWorkers = atoi(mammut::utils::split(lines.at(i), '\t')[0].c_str());
        if(currNumWorkers > numWorkers){
            numWorkers = currNumWorkers;
        }
    }

    // Farm is not used, we just create it to let things work.
    ff::ff_farm<> farm;
    DummyNode E;
    farm.add_emitter(&E);
    std::vector<ff::ff_node*> w;
    for(size_t i = 0;i < numWorkers; ++i) w.push_back(new DummyNode);
    farm.add_workers(w);
    DummyNode C;
    farm.add_collector(&C);

    nornir::Parameters ap("parameters.xml");
    nornir::ManagerFarm<> amf(&farm, ap);
    amf.simulate(std::string(argv[1]));
    amf.join();
    return 0;
}

