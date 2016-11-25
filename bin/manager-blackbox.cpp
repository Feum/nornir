/*
 * manager-external.cpp
 *
 * Created on: 21/06/2016
 *
 * This executable starts a manager which can monitor applications not written
 * with the Nornir framework.
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

#include <vector>
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "../src/manager.hpp"

using namespace nornir;

#undef DEBUG
#undef DEBUGB

#define DEBUG_MBB

#ifdef DEBUG_MBB
#define DEBUG(x) do { cerr << "[BlackBox Manager] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

static volatile bool *started, *handlerCreated;

int main(int argc, char * argv[]){
    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " ParametersFile Executable [Args]" << std::endl;
        return -1;
    }
    started = (volatile bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                                    MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    handlerCreated = (volatile bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                                           MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    *started = false;
    *handlerCreated = false;

    pid_t pid = fork();

    if(pid == -1){
        std::cerr << "Fork failed." << std::endl;
        return -1;
    }else if(pid){
        // Manager
        Parameters p(argv[1]);
        Observer o;
        p.observer = &o;
        std::cout << "Waiting for start. " << std::endl;
        while(!*started){;}
        std::cout << "Started." << std::endl;

        ManagerBlackBox m(p.mammut.getInstanceTask()->getProcessHandler(pid), p);
        *handlerCreated = true;
        m.start();
        m.join();
        waitpid(pid, NULL, 0);
    }else{
        // Application
        *started = true;
        // If the handler is not created, we would not catch the counters of
        // the threads/processes created by this process
        while(!*handlerCreated){;}
        munmap((void*) started, sizeof(volatile bool));
        munmap((void*) handlerCreated, sizeof(volatile bool));
        extern char** environ;
        if(execve(argv[2], &(argv[2]), environ) == -1){
            std::cerr << "Impossible to run the specified executable." << std::endl;
            return -1;
        }
        return -1; // execve never returns (except in the error case above).
    }
}

