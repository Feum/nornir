/*
 * utils.cpp
 *
 * Created on: 09/07/2015
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

/*!
 * \file utils.cpp
 * \brief Implementation of various utilities.
 **/

#include "utils.hpp"

namespace nornir{

using namespace std;
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::topology;

double Observer::calibrationDurationToPerc(const CalibrationStats& cs,
                                           uint durationMs){
    return ((double)cs.duration / (double)durationMs) * 100.0;
}

Observer::Observer(string statsFile, string calibrationFile, string summaryFile):
        _startMonitoringMs(0),
        _totalJoules(0),
        _lastTimestamp(0){
    _statsFile.open(statsFile.c_str());
    _calibrationFile.open(calibrationFile.c_str());
    _summaryFile.open(summaryFile.c_str());
    if(!_statsFile ||
       !_calibrationFile ||
       !_summaryFile){
        throw runtime_error("Observer: Impossible to open file.");
    }
    _statsFile << "TimestampMillisecs" << "\t";
    _statsFile << "[VirtualCores]" << "\t";
    _statsFile << "Workers" << "\t";
    _statsFile << "Frequency" << "\t";
    _statsFile << "CurrentBandwidth" << "\t";
    _statsFile << "SmoothedBandwidth" << "\t";
    _statsFile << "CoeffVarBandwidth" << "\t";
    _statsFile << "SmoothedLatency" << "\t";
    _statsFile << "SmoothedUtilization" << "\t";
    _statsFile << "CurrentWatts" << "\t";
    _statsFile << "SmoothedWatts" << "\t";
    _statsFile << endl;

    _calibrationFile << "NumSteps" << "\t";
    _calibrationFile << "TimeMs" << "\t";
    _calibrationFile << "Time%" << "\t";
    _calibrationFile << "TasksNum" << "\t";
    _calibrationFile << "Tasks%" << "\t";
    _calibrationFile << "Watts" << "\t";
    _calibrationFile << endl;

    _summaryFile << "Watts" << "\t";
    _summaryFile << "Bandwidth" << "\t";
    _summaryFile << "CompletionTimeSec" << "\t";
    _summaryFile << "CalibrationSteps" << "\t";
    _summaryFile << "CalibrationTimeMs" << "\t";
    _summaryFile << "CalibrationTime%" << "\t";
    _summaryFile << "CalibrationTasksNum" << "\t";
    _summaryFile << "CalibrationTasks%" << "\t";
    _summaryFile << "CalibrationWatts" << "\t";
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        _summaryFile << "Reconfigurations" << knobTypeToString((KnobType) i) << "Average" << "\t";
        _summaryFile << "Reconfigurations" << knobTypeToString((KnobType) i) << "Stddev" << "\t";
    }
    _summaryFile << "ReconfigurationsTotalAverage" << "\t";
    _summaryFile << "ReconfigurationsTotalStddev" << "\t";
    _summaryFile << endl;
}

Observer::~Observer(){
    _statsFile.close();
    _calibrationFile.close();
    _summaryFile.close();
}

void Observer::addJoules(Joules j){
    _totalJoules += j;
}

void Observer::observe(unsigned int timeStamp,
                     size_t workers,
                     Frequency frequency,
                     const vector<VirtualCore*>& virtualCores,
                     double currentBandwidth,
                     double smoothedBandwidth,
                     double coeffVarBandwidth,
                     double smoothedLatency,
                     double smoothedUtilization,
                     Joules currentWatts,
                     Joules smoothedWatts){
    if(_lastTimestamp == 0){
        _lastTimestamp = _startMonitoringMs;
    }
    _lastTimestamp = timeStamp;
    _statsFile << timeStamp - _startMonitoringMs << "\t";
    _statsFile << "[";
    for(size_t i = 0; i < virtualCores.size(); i++){
        _statsFile << virtualCores.at(i)->getVirtualCoreId() << ",";
    }
    _statsFile << "]" << "\t";

    _statsFile << workers << "\t";
    _statsFile << frequency << "\t";
    _statsFile << currentBandwidth << "\t";
    _statsFile << smoothedBandwidth << "\t";
    _statsFile << coeffVarBandwidth << "\t";
    _statsFile << smoothedLatency << "\t";
    _statsFile << smoothedUtilization << "\t";

    _statsFile << currentWatts << "\t";
    _statsFile << smoothedWatts << "\t";

    _statsFile << endl;
}

void Observer::calibrationStats(const vector<CalibrationStats>& calibrationStats,
                                uint durationMs,
                                uint64_t totalTasks){

    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        _calibrationFile << cs.numSteps << "\t";
        _calibrationFile << cs.duration << "\t";
        _calibrationFile << calibrationDurationToPerc(cs, durationMs) << "\t";
        _calibrationFile << cs.numTasks << "\t";
        _calibrationFile << ((double) cs.numTasks / totalTasks) * 100.0 << "\t";
        _calibrationFile << cs.joules / (cs.duration / 1000.0)<< "\t";
        _calibrationFile << endl;
    }
}

void Observer::summaryStats(const vector<CalibrationStats>& calibrationStats,
                            ReconfigurationStats reconfigurationStats,
                            uint durationMs,
                            uint64_t totalTasks){

    CalibrationStats totalCalibration;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        totalCalibration += calibrationStats.at(i);
    }
    
    if(durationMs == totalCalibration.duration){
        durationMs += 0.001; // Just to avoid division by 0 
    }

    _summaryFile << (_totalJoules - totalCalibration.joules) / (double) ((durationMs - totalCalibration.duration) / 1000.0) << "\t";
    _summaryFile << (totalTasks - totalCalibration.numTasks) / (double) ((durationMs - totalCalibration.duration) / 1000.0) << "\t";
    _summaryFile << (double) durationMs / 1000.0 << "\t";
    _summaryFile << totalCalibration.numSteps << "\t";
    _summaryFile << totalCalibration.duration << "\t";
    _summaryFile << calibrationDurationToPerc(totalCalibration, durationMs) << "\t";
    _summaryFile << totalCalibration.numTasks << "\t";
    _summaryFile << ((double) totalCalibration.numTasks / totalTasks) * 100.0 << "\t";
    _summaryFile << totalCalibration.joules / (totalCalibration.duration / 1000.0) << "\t";
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(reconfigurationStats.storedKnob((KnobType) i)){
            _summaryFile << reconfigurationStats.getAverageKnob((KnobType) i) << "\t";
            _summaryFile << reconfigurationStats.getStdDevKnob((KnobType) i) << "\t";
        }else{
            _summaryFile << "N.D." << "\t";
            _summaryFile << "N.D." << "\t";
        }
    }

    if(reconfigurationStats.storedTotal()){
        _summaryFile << reconfigurationStats.getAverageTotal() << "\t";
        _summaryFile << reconfigurationStats.getStdDevTotal() << "\t";
    }else{
        _summaryFile << "N.D." << "\t";
        _summaryFile << "N.D." << "\t";
    }


    _summaryFile << endl;
}

}
