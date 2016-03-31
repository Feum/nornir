/*
 * stream.cpp
 *
 * Created on: 26/03/2016
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

#include "stream.hpp"

#include "../external/Mammut/mammut/mammut.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_DF_STREAM
#define DEBUG(x) do { cerr << "[Dataflow:Stream] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{
namespace dataflow{

using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;
using namespace std;


ClockThread::ClockThread(time_t& lastSec, bool& terminated):
    _lastSec(lastSec), _terminated(terminated){
    ;
}

void ClockThread::run(){
    Mammut m;
    TasksManager* pm = m.getInstanceTask();

    ProcessHandler* process = pm->getProcessHandler(getpid());
    ThreadHandler* thread = process->getThreadHandler(gettid());
    thread->move((VirtualCoreId) 0);
    process->releaseThreadHandler(thread);
    pm->releaseProcessHandler(process);

    _lastSec = time(NULL);
    while(!_terminated){
        sleep(1);
        _lastSec = time(NULL);
    }
}

InputStreamRate::InputStreamRate(const std::string& fileName):
        _currentInterval(0),
        _lastSec(0), _startTime(0), _processedPkts(0), _currIntervalStart(0),
        _currBurstSize(0), _excess(0), _def(0), _terminated(false),
        _nextObject(0){
    _clockThread = new ClockThread(_lastSec, _terminated);
    FILE* f = NULL;
    char line[512];
    f = fopen(fileName.c_str(), "r");
    float rate = 0;
    float duration = 0;
    if(f){
        while(fgets(line, 512, f) != NULL){
            sscanf(line, "%f %f", &rate, &duration);
            Rates r;
            r.rate = rate;
            r.duration = duration;
            _rates.push_back(r);
        }
        fclose(f);
    }

    Mammut m;
    CpuFreq* frequency = m.getInstanceCpuFreq();
    vector<Domain*> domains = frequency->getDomains();
    _clockFrequency = 0;
    for(size_t i = 0; i < domains.size(); i++){
        Domain* domain = domains.at(i);
        double maxFrequency = domain->getAvailableFrequencies().back() * 1000; // We need the frequency in Hz (mammut returns frequency in kHz)
        if(!_clockFrequency || maxFrequency == _clockFrequency){
            _clockFrequency = maxFrequency;
        }else{
            throw std::runtime_error("ERROR: Different max frequencies for different domains.");
        }
    }
}

InputStreamRate::~InputStreamRate(){
    delete _clockThread;
}

void InputStreamRate::init(){
    _objects = loadObjects();
    DEBUG("Loaded " << _objects.size() << " objects.");
    DEBUG(_rates.size() << " rates.");
}

#define BURST_SIZE 10.0

StreamElem* InputStreamRate::next(){
    if(!_objects.size()){
        throw std::runtime_error("You need to call init() before calling next() for the first time.");
    }

    StreamElem* obj = NULL;
    if(!_def){
        _def = getticks();
        _startTime = time(NULL);
        _clockThread->start();
    }

    if(_nextObject == _objects.size()){
        _nextObject = 0;
    }

    if(_currentInterval < _rates.size()){
        obj = _objects.at(_nextObject);
    }else{
        DEBUG("Stream terminated.");
        obj = NULL;
        _terminated = true;
        _clockThread->join();
        return obj;
    }

    if(_currBurstSize == BURST_SIZE){
        /** Sleep to get the rate. **/
        double wait_interval_secs = 1.0 / _rates[_currentInterval].rate;
        ticks ticks_to_sleep = (_clockFrequency * wait_interval_secs * (double) BURST_SIZE);

        DEBUG("Rate: " << _rates[_currentInterval].rate << " Wait interval secs: " << wait_interval_secs << " TTS: " << ticks_to_sleep);

        _currBurstSize = 0;

        _excess += (getticks()-_def);

        if(_excess >= ticks_to_sleep){
            //_excess = 0;
            _excess -= ticks_to_sleep;
        }else{
            _excess = ticksWait(ticks_to_sleep - _excess);
        }

        _def = getticks();
    }


    ++_currBurstSize;

    if(_currIntervalStart == 0){
        _currIntervalStart = _lastSec;
    }

    ++_processedPkts;

    /** Go to the next rate **/
    if(_lastSec - _currIntervalStart >= _rates[_currentInterval].duration){
        _currIntervalStart = _lastSec;
        _currentInterval++;
        DEBUG("Moving to interval: " <<  _currentInterval);
        _excess = 0;
    }

    ++_nextObject;

    return obj;
}

bool InputStreamRate::hasNext(){
    return _currentInterval < _rates.size();
}

}
}
