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
    _calibrationFile << "Duration" << "\t";
    _calibrationFile << "Time%" << "\t";
    _calibrationFile << endl;

    _summaryFile << "Watts" << "\t";
    _summaryFile << "Bandwidth" << "\t";
    _summaryFile << "CompletionTime" << "\t";
    _summaryFile << "Calibration%" << "\t";
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

void Observer::calibrationStats(const vector<CalibrationStats>&
                              calibrationStats,
                              uint durationMs){

    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        _calibrationFile << cs.numSteps << "\t";
        _calibrationFile << cs.duration << "\t";
        _calibrationFile << calibrationDurationToPerc(cs, durationMs) << "\t";
        _calibrationFile << endl;
    }
}

void Observer::summaryStats(const vector<CalibrationStats>&
                          calibrationStats,
                          uint durationMs){
    double totalCalibrationPerc = 0.0;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        totalCalibrationPerc += calibrationDurationToPerc(cs, durationMs);
    }

    _summaryFile << _totalJoules / (double) (durationMs / 1000.0) << "\t";
    _summaryFile << _totalBw / (double) _numSamples << "\t";
    _summaryFile << (double) durationMs / 1000.0 << "\t";
    _summaryFile << totalCalibrationPerc << "\t";
    _summaryFile << endl;
}

}
