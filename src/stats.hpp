/*
 * stats.hpp
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
 * \file stats.hpp
 * \brief Statistics collection utilities.
 **/

#ifndef NORNIR_STATS_HPP_
#define NORNIR_STATS_HPP_

#include "knob.hpp"
#include "utils.hpp"

namespace nornir{


class ReconfigurationStats{
private:
    std::vector<double> _knobs[KNOB_TYPE_NUM];
    std::vector<double> _total;
    bool _storedKnob[KNOB_TYPE_NUM];
    bool _storedTotal;
public:
    ReconfigurationStats(){
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            _storedKnob[i] = false;
        }
        _storedTotal = false;
    }

    void swap(ReconfigurationStats& x){
        using std::swap;
        swap(_knobs, x._knobs);
        swap(_total, x._total);
        swap(_storedKnob, x._storedKnob);
        swap(_storedTotal, x._storedTotal);
    }

    inline ReconfigurationStats(const ReconfigurationStats& other){
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            _knobs[i] = other._knobs[i];
            _storedKnob[i] = other._storedKnob[i];
        }
        _total = other._total;
        _storedTotal = other._storedTotal;
    }

    inline ReconfigurationStats& operator=(ReconfigurationStats other){
        swap(other);
        return *this;
    }

    inline void addSample(KnobType idx, double sample){
        _storedKnob[idx] = true;
        _knobs[idx].push_back(sample);
    }

    inline void addSampleTotal(double total){
        _storedTotal = true;
        _total.push_back(total);
    }

    inline double getAverageKnob(KnobType idx){
        return average(_knobs[idx]);
    }

    inline double getStdDevKnob(KnobType idx){
        return stddev(_knobs[idx]);
    }

    inline double getAverageTotal(){
        return average(_total);
    }

    inline double getStdDevTotal(){
        return stddev(_total);
    }

    inline bool storedKnob(KnobType idx){
        return _storedKnob[idx];
    }

    inline bool storedTotal(){
        return _storedTotal;
    }
};

typedef struct CalibrationStats{
    uint numSteps;
    uint duration;
    uint numTasks;
    mammut::energy::Joules joules;

    CalibrationStats():numSteps(0), duration(0),numTasks(0),joules(0){
        ;
    }

    void swap(CalibrationStats& x){
        using std::swap;

        swap(numSteps, x.numSteps);
        swap(duration, x.duration);
        swap(numTasks, x.numTasks);
        swap(joules, x.joules);
    }

    CalibrationStats& operator=(CalibrationStats rhs){
        swap(rhs);
        return *this;
    }

    CalibrationStats& operator+=(const CalibrationStats& rhs){
        numSteps += rhs.numSteps;
        duration += rhs.duration;
        numTasks += rhs.numTasks;
        joules += rhs.joules;
        return *this;
    }
}CalibrationStats;

/*!
 * This class can be used to obtain statistics about reconfigurations
 * performed by the manager.
 * It can be extended by a user defined class to customize action to take
 * for each observed statistic.
 */
class Observer{
    friend class Manager;
private:
    std::ofstream _statsFile;
    std::ofstream _calibrationFile;
    std::ofstream _summaryFile;
    unsigned int _startMonitoringMs;
    mammut::energy::Joules _totalJoules;
    unsigned int _lastTimestamp;

    void addJoules(mammut::energy::Joules j);

    double calibrationDurationToPerc(const CalibrationStats& cs,
                                     uint durationMs);
public:
    Observer(std::string statsFile = "stats.csv",
             std::string calibrationFile = "calibration.csv",
             std::string summaryFile = "summary.csv");

    virtual ~Observer();

    virtual void observe(unsigned int timeStamp,
                         size_t workers,
                         mammut::cpufreq::Frequency frequency,
                         const std::vector<mammut::topology::VirtualCore*>& virtualCores,
                         double currentBandwidth,
                         double smoothedBandwidth,
                         double coeffVarBandwidth,
                         double smoothedLatency,
                         double smoothedUtilization,
                         mammut::energy::Joules currentWatts,
                         mammut::energy::Joules smoothedWatts);

    virtual void calibrationStats(const std::vector<CalibrationStats>&
                                  calibrationStats,
                                  uint durationMs,
                                  uint64_t totalTasks);

    virtual void summaryStats(const std::vector<CalibrationStats>& calibrationStats,
                              ReconfigurationStats reconfigurationStats,
                              uint durationMs,
                              uint64_t totalTasks);
};

}

#endif
