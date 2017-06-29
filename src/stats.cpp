/*
 * stats.cpp
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
 * \file stats.cpp
 * \brief Statistics collection utilities.
 **/

#include "stats.hpp"
#include "configuration.hpp"
#include "selectors.hpp"
extern "C"{
#include "external/graphite-c-client/graphite-client.h"
}

namespace nornir{

using namespace std;
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::topology;

Logger::Logger(unsigned int timeOffset):
    _timeOffset(timeOffset), _startMonitoring(0), _totalJoules(0){;}

void Logger::addJoules(Joules j){
    _totalJoules += j;
}

unsigned int Logger::getAbsoluteTimestamp(){
    return mammut::utils::getMillisecondsTime();
}

unsigned int Logger::getRelativeTimestamp(){
    unsigned int absolute = getAbsoluteTimestamp();
    if(!_startMonitoring){
        _startMonitoring = absolute;
    }

    return absolute - _startMonitoring + _timeOffset;
}

LoggerStream::LoggerStream(std::ostream *statsStream,
                           std::ostream *calibrationStream,
                           std::ostream *summaryStream,
                           unsigned int timeOffset):
        Logger(timeOffset),
        _statsStream(statsStream), _calibrationStream(calibrationStream),
        _summaryStream(summaryStream), _timeOffset(timeOffset){
    if(!*_statsStream ||
       !*_calibrationStream ||
       !*_summaryStream){
        throw runtime_error("LoggerOutStream: Impossible to use stream.");
    }
    *_statsStream << "TimestampMillisecs" << "\t";
    *_statsStream << "[VirtualCores]" << "\t";
    *_statsStream << "Workers" << "\t";
    *_statsStream << "Frequency" << "\t";
    *_statsStream << "CurrentBandwidth" << "\t";
    *_statsStream << "SmoothedBandwidth" << "\t";
    *_statsStream << "CoeffVarBandwidth" << "\t";
    *_statsStream << "SmoothedLatency" << "\t";
    *_statsStream << "SmoothedUtilization" << "\t";
    *_statsStream << "CurrentWatts" << "\t";
    *_statsStream << "SmoothedWatts" << "\t";
    *_statsStream << endl;

    *_calibrationStream << "NumSteps" << "\t";
    *_calibrationStream << "TimeMs" << "\t";
    *_calibrationStream << "Time%" << "\t";
    *_calibrationStream << "TasksNum" << "\t";
    *_calibrationStream << "Tasks%" << "\t";
    *_calibrationStream << "Watts" << "\t";
    *_calibrationStream << endl;

    *_summaryStream << "Watts" << "\t";
    *_summaryStream << "Bandwidth" << "\t";
    *_summaryStream << "CompletionTimeSec" << "\t";
    *_summaryStream << "CalibrationSteps" << "\t";
    *_summaryStream << "CalibrationTimeMs" << "\t";
    *_summaryStream << "CalibrationTime%" << "\t";
    *_summaryStream << "CalibrationTasksNum" << "\t";
    *_summaryStream << "CalibrationTasks%" << "\t";
    *_summaryStream << "CalibrationWatts" << "\t";
    for(size_t i = 0; i < KNOB_NUM; i++){
        *_summaryStream << "Reconfigurations" << knobTypeToString((KnobType) i) << "Average" << "\t";
        *_summaryStream << "Reconfigurations" << knobTypeToString((KnobType) i) << "Stddev" << "\t";
    }
    *_summaryStream << "ReconfigurationsTotalAverage" << "\t";
    *_summaryStream << "ReconfigurationsTotalStddev" << "\t";
    *_summaryStream << endl;
}

void LoggerStream::log(const Configuration& configuration, const Smoother<MonitoredSample>& samples){
    const vector<VirtualCore*>& virtualCores = dynamic_cast<const KnobMapping*>(configuration.getKnob(KNOB_MAPPING))->getActiveVirtualCores();
    MonitoredSample ms = samples.average();
    *_statsStream << getRelativeTimestamp() << "\t";
    *_statsStream << "[";
    for(size_t i = 0; i < virtualCores.size(); i++){
        *_statsStream << virtualCores.at(i)->getVirtualCoreId() << ",";
    }
    *_statsStream << "]" << "\t";

    *_statsStream << configuration.getRealValue(KNOB_VIRTUAL_CORES) << "\t";
    // Print frequency as string to avoid conversion to exp notation.
    std::ostringstream strs;
    strs << std::fixed << std::setprecision(0) << configuration.getRealValue(KNOB_FREQUENCY);
    *_statsStream << strs.str() << "\t";
    *_statsStream << samples.getLastSample().bandwidth << "\t";
    *_statsStream << ms.bandwidth << "\t";
    *_statsStream << samples.coefficientVariation().bandwidth << "\t";
    *_statsStream << ms.latency << "\t";
    *_statsStream << ms.loadPercentage << "\t";

    *_statsStream << samples.getLastSample().watts << "\t";
    *_statsStream << ms.watts << "\t";

    *_statsStream << endl;
}

void LoggerStream::logSummary(const Configuration& configuration, Selector* selector, ulong durationMs, double totalTasks){
    vector<CalibrationStats> calibrationStats;
    ReconfigurationStats reconfigurationStats = configuration.getReconfigurationStats();
    if(selector){
        selector->updateTotalTasks(totalTasks);
        selector->stopCalibration();
        calibrationStats = selector->getCalibrationsStats();
    }

    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        *_calibrationStream << cs.numSteps << "\t";
        *_calibrationStream << cs.duration << "\t";
        *_calibrationStream << ((double) cs.duration / (double)durationMs) * 100.0 << "\t";
        *_calibrationStream << cs.numTasks << "\t";
        *_calibrationStream << ((double) cs.numTasks / totalTasks) * 100.0 << "\t";
        *_calibrationStream << cs.joules / (cs.duration / 1000.0)<< "\t";
        *_calibrationStream << endl;
    }

    CalibrationStats totalCalibration;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        totalCalibration += calibrationStats.at(i);
    }

    if(durationMs == totalCalibration.duration){
        durationMs += 0.001; // Just to avoid division by 0
    }

    *_summaryStream << (_totalJoules - totalCalibration.joules) / (double) ((durationMs - totalCalibration.duration) / 1000.0) << "\t";
    *_summaryStream << (totalTasks - totalCalibration.numTasks) / (double) ((durationMs - totalCalibration.duration) / 1000.0) << "\t";
    *_summaryStream << (double) durationMs / 1000.0 << "\t";
    *_summaryStream << totalCalibration.numSteps << "\t";
    *_summaryStream << totalCalibration.duration << "\t";
    *_summaryStream << ((double) totalCalibration.duration / (double)durationMs) * 100.0 << "\t";
    *_summaryStream << totalCalibration.numTasks << "\t";
    *_summaryStream << ((double) totalCalibration.numTasks / totalTasks) * 100.0 << "\t";
    *_summaryStream << totalCalibration.joules / (totalCalibration.duration / 1000.0) << "\t";
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(reconfigurationStats.storedKnob((KnobType) i)){
            *_summaryStream << reconfigurationStats.getAverageKnob((KnobType) i) << "\t";
            *_summaryStream << reconfigurationStats.getStdDevKnob((KnobType) i) << "\t";
        }else{
            *_summaryStream << "N.D." << "\t";
            *_summaryStream << "N.D." << "\t";
        }
    }

    if(reconfigurationStats.storedTotal()){
        *_summaryStream << reconfigurationStats.getAverageTotal() << "\t";
        *_summaryStream << reconfigurationStats.getStdDevTotal() << "\t";
    }else{
        *_summaryStream << "N.D." << "\t";
        *_summaryStream << "N.D." << "\t";
    }

    *_summaryStream << endl;
}

LoggerGraphite::LoggerGraphite(const std::string& host, unsigned int port){
    graphite_init(host.c_str(), port);
}

LoggerGraphite::~LoggerGraphite(){
    graphite_finalize();
}

void LoggerGraphite::log(const Configuration& configuration,
                       const Smoother<MonitoredSample>& samples){
    unsigned int timestamp = getAbsoluteTimestamp();


    /*************************************************/
    /*               Resources info                  */
    /*************************************************/
    for(auto vc : dynamic_cast<const KnobMapping*>(configuration.getKnob(KNOB_MAPPING))->getActiveVirtualCores()){
        graphite_send_plain((std::string("nornir.resources.cores") + utils::intToString(vc->getVirtualCoreId())).c_str(), 1, timestamp);
    }
    // Set all other cores to zero.
    for(auto vc : dynamic_cast<const KnobMapping*>(configuration.getKnob(KNOB_MAPPING))->getUnusedVirtualCores()){
        graphite_send_plain((std::string("nornir.resources.cores") + utils::intToString(vc->getVirtualCoreId())).c_str(), 0, timestamp);
    }
    graphite_send_plain("nornir.resources.cores.num", configuration.getRealValue(KNOB_VIRTUAL_CORES), timestamp);
    graphite_send_plain("nornir.resources.frequency", configuration.getRealValue(KNOB_FREQUENCY), timestamp);

    /*************************************************/
    /*                Monitor info                  */
    /*************************************************/
    graphite_send_plain("nornir.monitor.bandwidth.current", samples.getLastSample().bandwidth, timestamp);
    graphite_send_plain("nornir.monitor.bandwidth.average", samples.average().bandwidth, timestamp);
    graphite_send_plain("nornir.monitor.power.current", samples.getLastSample().watts, timestamp);
    graphite_send_plain("nornir.monitor.power.average", samples.average().watts, timestamp);
    graphite_send_plain("nornir.monitor.latency.current", samples.getLastSample().latency, timestamp);
    graphite_send_plain("nornir.monitor.latency.average", samples.average().latency, timestamp);
    graphite_send_plain("nornir.monitor.utilization.current", samples.getLastSample().loadPercentage, timestamp);
    graphite_send_plain("nornir.monitor.utilization.average", samples.average().loadPercentage, timestamp);
}

}
