/*
 * farm_ffcompat.hpp
 *
 * Created on: 27/02/2016
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
 * Basic test for the nornir farm.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "../src/interface.hpp"

using namespace ff;

#define MICROSECSSLEEP 1000000

static int maxTasks;

/**
 * Scheduler.
 */
class Emitter: public nornir::Scheduler<int>{
public:
    int* schedule() {
        usleep(MICROSECSSLEEP);
        int * task = new int(maxTasks);
        --maxTasks;
        if (maxTasks < 0){
            std::cout << "Scheduler finished" << std::endl;
            return NULL;
        }
        return task;
    }
};

/**
 * Worker.
 */
class Worker: public nornir::Worker<int, int>{
public:
    int * compute(int * task) {
        usleep(MICROSECSSLEEP);
        std::cout << "Worker " << getId()
                  << " received task " << *task << std::endl;
        return task;
    }
};

/**
 * Gatherer.
 */
class Collector: public nornir::Gatherer<int> {
public:
    void gather(int* task) {
        std::cout << "Gatherer received task " << *task << std::endl;
    }
};

int main(int argc, char * argv[]) {
    int nworkers = 1;
    int streamlen = 10;

    if (argc>1) {
        if (argc!=3) {
            std::cerr << "use: " 
                      << argv[0] 
                      << " nworkers streamlen\n";
            return -1;
        }   
        nworkers = atoi(argv[1]);
        streamlen = atoi(argv[2]);
    }

    if (!nworkers || !streamlen) {
        std::cerr << "Wrong parameters values\n";
        return -1;
    }

    maxTasks = streamlen;

    nornir::Farm<int, int> farm("parameters.xml");
    farm.start<Emitter, Worker, Collector>(nworkers);
    farm.wait();

    return 0;
}

