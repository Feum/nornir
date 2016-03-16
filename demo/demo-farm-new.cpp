/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  As a special exception, you may use this file as part of a free software
 *  library without restriction.  Specifically, if other files instantiate
 *  templates or use macros or inline functions from this file, or you compile
 *  this file and link it with other files to produce an executable, this
 *  file does not by itself cause the resulting executable to be covered by
 *  the GNU General Public License.  This exception does not however
 *  invalidate any other reasons why the executable file might be covered by
 *  the GNU General Public License.
 *
 ****************************************************************************
 */

/**

 Very basic test for the FastFlow farm.
 
*/
#include <vector>
#include <iostream>
#include <fstream>
#include "../src/interface.hpp"

using namespace ff;

//#define MICROSECSSLEEP 200000
#define MICROSECSSLEEP 1000000

static int maxTasks;

class Emitter: public nornir::Scheduler<int>{
public:
    int* schedule() {
        usleep(MICROSECSSLEEP);
        int * task = new int(maxTasks);
        --maxTasks;
        if (maxTasks < 0){
            std::cout << "Emitter finished" << std::endl;
            return NULL;
        }
        return task;
    }
};


// generic worker
class Worker: public nornir::Worker<int, int>{
public:
    int * compute(int * task) {
        usleep(MICROSECSSLEEP);
        std::cout << "Worker " << ff_node::get_my_id()
                  << " received task " << *task << "\n";
        return task;
    }
};

// the gatherer filter
class Collector: public nornir::Gatherer<int> {
public:
    void gather(int* task) {
        std::cout << "Collector received task " << *task << "\n";
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

    nornir::Farm<int, int> farm("parameters.xml", "archdata.xml");
    std::cout << "Starting farm. " << std::endl;
    farm.start<Emitter, Worker, Collector>(nworkers);
    std::cout << "Farm started. " << std::endl;
    farm.wait();
    std::cout << "Farm joined. " << std::endl;

    return 0;
}

