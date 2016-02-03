/*
 * utils.cpp
 *
 * Created on: 09/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of AdaptiveFastFlow.
 *
 *  AdaptiveFastFlow is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  AdaptiveFastFlow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with AdaptiveFastFlow.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

/*!
 * \file utils.cpp
 * \brief Implementation of various utilities.
 **/

#include "utils.hpp"
#include <cmath>
#include <functional>
#include <numeric>

namespace adpff{

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
        _totalBw(0),
        _numSamples(0),
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
    _statsFile << "[[EmitterVc][WorkersVc][CollectorVc]]" << "\t";
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
    _calibrationFile << endl;

    _summaryFile << "Watts" << "\t";
    _summaryFile << "Bandwidth" << "\t";
    _summaryFile << "CompletionTimeSec" << "\t";
    _summaryFile << "CalibrationSteps" << "\t";
    _summaryFile << "CalibrationTimeMs" << "\t";
    _summaryFile << "CalibrationTime%" << "\t";
    _summaryFile << "CalibrationTasksNum" << "\t";
    _summaryFile << "CalibrationTasks%" << "\t";
    _summaryFile << "ReconfigurationsMsAverage" << "\t";
    _summaryFile << "ReconfigurationsMsStddev" << "\t";
    _summaryFile << endl;
}

Observer::~Observer(){
    _statsFile.close();
    _calibrationFile.close();
    _summaryFile.close();
}

void Observer::observe(unsigned int timeStamp,
                     size_t workers,
                     Frequency frequency,
                     const VirtualCore* emitterVirtualCore,
                     const vector<VirtualCore*>& workersVirtualCore,
                     const VirtualCore* collectorVirtualCore,
                     double currentBandwidth,
                     double smoothedBandwidth,
                     double coeffVarBandwidth,
                     double smoothedLatency,
                     double smoothedUtilization,
                     Joules currentWatts,
                     Joules smoothedWatts){
    unsigned int interval;
    if(_lastTimestamp == 0){
        _lastTimestamp = _startMonitoringMs;
    }
    interval = timeStamp - _lastTimestamp;
    _lastTimestamp = timeStamp;
    _statsFile << timeStamp - _startMonitoringMs << "\t";
    _statsFile << "[";
    if(emitterVirtualCore){
        _statsFile << "[" << emitterVirtualCore->getVirtualCoreId() << "]";
    }

    _statsFile << "[";
    for(size_t i = 0; i < workersVirtualCore.size(); i++){
        _statsFile << workersVirtualCore.at(i)->getVirtualCoreId() << ",";
    }
    _statsFile << "]";

    if(collectorVirtualCore){
        _statsFile << "[" << collectorVirtualCore->getVirtualCoreId() << "]";
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

    _totalJoules += currentWatts * (interval / 1000.0);
    _totalBw += currentBandwidth;
    _numSamples++;
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
        _calibrationFile << endl;
    }
}

inline double average(const vector<double>& v){
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

inline double stddev(const vector<double>& v, double average){
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   bind2nd(std::minus<double>(), average));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return sqrt(sq_sum / v.size());
}

inline double stddev(const vector<double>& v){
    return stddev(v, average(v));
}


void Observer::summaryStats(const vector<CalibrationStats>& calibrationStats,
                            const vector<double>& reconfigurationStats,
                            uint durationMs,
                            uint64_t totalTasks){

    CalibrationStats totalCalibration;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        totalCalibration += calibrationStats.at(i);
    }

    _summaryFile << _totalJoules / (double) (durationMs / 1000.0) << "\t";
    _summaryFile << _totalBw / (double) _numSamples << "\t";
    _summaryFile << (double) durationMs / 1000.0 << "\t";
    _summaryFile << totalCalibration.numSteps << "\t";
    _summaryFile << totalCalibration.duration << "\t";
    _summaryFile << calibrationDurationToPerc(totalCalibration, durationMs) << "\t";
    _summaryFile << totalCalibration.numTasks << "\t";
    _summaryFile << ((double) totalCalibration.numTasks / totalTasks) * 100.0 << "\t";
    if(reconfigurationStats.size()){
        double reconfigurationMsAvg = average(reconfigurationStats);
        _summaryFile << reconfigurationMsAvg << "\t";
        _summaryFile << stddev(reconfigurationStats, reconfigurationMsAvg) << "\t";
    }else{
        _summaryFile << "N.D." << "\t";
        _summaryFile << "N.D." << "\t";
    }
    _summaryFile << endl;
}

}
