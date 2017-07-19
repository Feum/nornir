/*
 * manager-external.cpp
 *
 * Created on: 21/06/2016
 *
 * This executable starts a manager which can monitor applications not written
 * with the Nornir framework. Albeit more of such applications can run concurrently,
 * no coordination between them is provided. The result may thus be inconsistent.
 * Please use this binary to control at most one application at a time.
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
#include <stdlib.h>
#include "../src/manager.hpp"
#include "../src/manager-multi.hpp"
#include "../src/monitor.hpp"

using namespace nornir;
using namespace nn;

#undef DEBUG
#undef DEBUGB

#define DEBUG_MEXT

#ifdef DEBUG_MEXT
#define DEBUG(x) do { cerr << "[External Server] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

class ApplicationInstance{
public:
    nn::socket channel;
    int chid;
    ManagerInstrumented *manager;

    ApplicationInstance():channel(AF_SP, NN_PAIR), chid(0), manager(NULL), observer(NULL){;}
};

int main(int argc, char * argv[]){
    nn::socket mainChannel(AF_SP, NN_PAIR);
    int mainChid;
    mainChid = mainChannel.bind(INSTRUMENTATION_CONNECTION_CHANNEL);
    std::list<ApplicationInstance*> instances;
    if(system("mkdir -p /tmp/nornir") == -1){throw std::runtime_error("Impossible to create nornir dir.");}
    if(system("chmod ugo+rwx /tmp/nornir.ipc") == -1){throw std::runtime_error("Impossible to set permission to nornir channel.");}
    if(system("chmod ugo+rwx /tmp/nornir/") == -1){throw std::runtime_error("Impossible to set permission nornir dir.");}
    ManagerMulti mm;
    mm.start();
    while(true){
        pid_t pid;
        size_t r = mainChannel.recv(&pid, sizeof(pid), 0);
        assert(r == sizeof(pid));

        DEBUG("Received a request from process " << pid);
        ApplicationInstance* ai = new ApplicationInstance;
        ai->chid = ai->channel.bind((string("ipc:///tmp/nornir/") + mammut::utils::intToString(pid) + string(".ipc")).c_str());
        DEBUG("Created app channel.");

        DEBUG("Sending ack.");
        char ack = 0;
        r = mainChannel.send(&ack, sizeof(ack), 0);
        assert(r == sizeof(ack));

        size_t length = 0;
        r = ai->channel.recv(&length, sizeof(length), 0);
        assert(r == sizeof(length));
        char* parameters = new char[length];
        DEBUG("Receiving parameters.");
        r = ai->channel.recv(parameters, length*sizeof(char), 0);
        assert(r == (int) (sizeof(char)*length));
        std::string parametersString(parameters);
        std::ofstream out("parameters.xml");
        out << parametersString;
        out.close();
        delete[] parameters;

        //TODO: If we will let this work for remote machines too, we will need
        // to also send archfile.xml
        Parameters p("parameters.xml");
        ConfigurationValidation pv = p.validate();
        DEBUG("Sending validation result.");
        r = ai->channel.send(&pv, sizeof(pv), 0);
        assert(r == sizeof(pv));
        ai->manager = new ManagerInstrumented(ai->channel, ai->chid, p);
        ai->manager->start();
        instances.push_back(ai);
        mm.addManager(ai->manager);

        /** Try to join and delete already terminated managers. **/
        Manager* m;
        m = mm.getTerminatedManager();
        if(m){
            for(auto it = instances.begin(); it != instances.end(); it++){
                if((*it)->manager == m){
                    auto newit = it + 1;
                    DEBUG("Application manager terminated, cleaning.");
                    delete ((*it)->manager);
                    (*it)->channel.shutdown((*it)->chid);
                    ApplicationInstance* ai = (*it);
                    instances.erase(it);
                    it = newit;
                    delete ai;
                }
            }
        }
    }
    mainChannel.shutdown(mainChid);
    return 0;
}

