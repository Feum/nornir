/*
 * monitor.cpp
 *
 * Created on: 03/07/2016
 *
 * This file contains the client code for monitoring a local application and sending
 * the monitoring data to an external manager (server).
 * Firs of all, we try to open the EXTERNAL_CHANNEL_NAME channel, whose name
 * is fixed and known by all the clients and by the server. On this channel,
 * the client sends its PID. Then, all the communications between this specific
 * client (i.e. application) and the server will be done on a separate channel
 * (identified by the PID).
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

#include <sys/types.h>
#include <unistd.h>

#include "parameters.hpp"
#include "monitor.hpp"
#include "external/mammut/mammut/mammut.hpp"

namespace nornir{

using namespace std;
using namespace mammut;
using namespace mammut::utils;

ExternalApplication::ExternalApplication(const std::string& parametersFile):
        _channel(AF_SP, NN_PAIR){

    /** Send pid, then switch to the pid channel. */
    nn::socket mainChannel(AF_SP, NN_PAIR);
    int mainChid;
    mainChid = mainChannel.connect(EXTERNAL_CHANNEL_NAME);
    if(mainChid < 0){
        throw std::runtime_error("Impossible to connect to external Nornir channel.");
    }
    pid_t pid = getpid();
    assert(mainChannel.send(&pid, sizeof(pid), 0) == sizeof(pid));
    // Wait ack. This is needed because we have to be sure that the pid channel
    // has been created before trying to connect.
    char ack;
    assert(mainChannel.recv(&ack, sizeof(ack), 0) == sizeof(ack));
    mainChannel.shutdown(mainChid);

    /** Send content length. **/
    std::string s;
    std::stringstream out;
    out << pid;
    _chid = _channel.connect((string("ipc:///tmp/nornir/") + out.str() + string(".ipc")).c_str());
    assert(_chid >= 0);
    vector<string> lines = readFile(parametersFile);
    size_t length = 0;
    for(size_t i = 0; i < lines.size(); i++){
        length += lines.at(i).length();
    }
    length += 1;
    assert(_channel.send(&length, sizeof(length), 0) == sizeof(length));

    /** Send content. **/
    string fileContent;
    fileContent.reserve(length);
    for(size_t i = 0; i < lines.size(); i++){
        fileContent.append(lines.at(i));
    }
    fileContent.push_back('\0');
    const char* fileContentC = fileContent.c_str();
    assert(_channel.send(fileContentC, length, 0) == (int) length);

    /** Receive validation result. **/
    ParametersValidation pv;
    assert(_channel.recv(&pv, sizeof(pv), 0) == sizeof(pv));
    if(pv != VALIDATION_OK){
        throw runtime_error("Invalid adaptivity parameters: " + std::to_string(pv));
    }

    /** Create orlog app. **/
    _app = new orlog::Application(_channel, _chid);
}

ExternalApplication::ExternalApplication(const std::string& parametersFile,
               const std::string& serverAddress,
               uint port):_channel(AF_SP, NN_PAIR), _chid(0), _app(NULL){
    throw std::runtime_error("Not yet implemented.");
}

ExternalApplication::~ExternalApplication(){
    _channel.shutdown(_chid);
    delete _app;
}

void ExternalApplication::begin(){
    _app->begin();
}

void ExternalApplication::end(){
    _app->end();
}

void ExternalApplication::terminate(){
    _app->terminate();
}

}
