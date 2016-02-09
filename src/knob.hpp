/*
 * knob.hpp
 *
 * Created on: 02/11/2015
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

#ifndef SRC_KNOB_HPP_
#define SRC_KNOB_HPP_

#include "node.hpp"
#include <ff/farm.hpp>

namespace adpff{

typedef enum{
    KNOB_TYPE_WORKERS = 0,
    KNOB_TYPE_MAPPING,
    KNOB_TYPE_FREQUENCY,
    KNOB_TYPE_NUM // <---- This must always be the last value
}KnobType;

class Knob: public mammut::utils::NonCopyable{
public:
    Knob():_realValue(-1){;}

    /**
     * Computes the real value corresponding to a specific
     * relative value.
     * @param relative The relative value.
     * @param real The real value.
     * @return True if there is a real value, false otherwise (e.g. because this
     * knob has no possible values to be set).
     */
    bool getRealFromRelative(double relative, double& real) const;

    /**
     * Changes the value of this knob.
     * @param v Is a value in the range [0, 100], representing the value to be
     *          set for the knob.
     */
    void setRelativeValue(double v);

    /**
     * Changes the value of this knob.
     * @param v Is the real value of the knob.
     */
    void setRealValue(double v);

    /**
     * Sets this knob to its maximum.
     */
    void setToMax();

    /**
     * Returns the current real value of this knob.
     * @return The current real value of this knob.
     */
    double getRealValue() const;

    /**
     * Returns true if this knob needs to be calibrated.
     * @return True if this knob needs to be calibrated.
     */
    virtual bool needsCalibration() const = 0;

    /**
     * Returns a vector of allowed values for this knob.
     * @return A vector of allowed values for this knob.
     */
    virtual std::vector<double> getAllowedValues() const = 0;

    virtual ~Knob(){;}
protected:
    /**
     * Changes the value of this knob.
     * @param v Is the real value that this knob will have.
     */
    virtual void changeValueReal(double v) = 0;

    double _realValue;
};

class KnobWorkers: public Knob{
public:
    KnobWorkers(KnobConfWorkers confWorkers, AdaptiveNode* emitter,
                AdaptiveNode* collector, ff::ff_gatherer* gt,
                const std::vector<AdaptiveNode*> workers,
                const volatile bool* terminated);

    bool needsCalibration() const;
    void changeValueReal(double v);
    std::vector<double> getAllowedValues() const;

    /**
     * Returns the number of cores on which the workers should be
     * executed.
     * @return The number of cores on which the workers should be
     * executed.
     */
    uint getNumActiveCores() const;

    /**
     * Returns a vector containing all the active workers.
     * @return A vector containing all the active workers.
     */
    const std::vector<AdaptiveNode*>& getActiveWorkers() const;
private:
    /**
     * Prepares the nodes to freeze.
     */
    void prepareToFreeze();

    /**
     * Freezes all the nodes.
     */
    void freeze();

    /**
     * Prepares the nodes to run.
     * @param numWorkers The new number of workers.
     */
    void prepareToRun(uint numWorkers);

    /**
     * Runs the farm with a new number of workers.
     * @param numWorkers The new number of workers.
     */
    void run(uint numWorkers);

    /**
     * Notifies a change in the number of workers to all the nodes.
     * @param numWorkers The new number of workers.
     */
    void notifyNewConfiguration(uint numWorkers);

    KnobConfWorkers _confWorkers;
    AdaptiveNode* _emitter;
    AdaptiveNode* _collector;
    ff::ff_gatherer* _gt;
    const std::vector<AdaptiveNode*> _allWorkers;
    std::vector<AdaptiveNode*> _activeWorkers;
    std::vector<double> _knobValues;
    const volatile bool* _terminated;
};

class KnobMapping: public Knob{
public:
    KnobMapping(KnobConfMapping confMapping,
                KnobConfSNodeMapping confEmitterMapping,
                KnobConfSNodeMapping confCollectorMapping,
                KnobConfHyperthreading confHyperthreading,
                const mammut::Mammut& mammut,
                AdaptiveNode* emitter,
                AdaptiveNode* collector,
                const KnobWorkers& knobWorkers);
    bool needsCalibration() const;
    void changeValueReal(double v);
    std::vector<double> getAllowedValues() const;

    mammut::topology::VirtualCore* getEmitterVirtualCore() const;
    mammut::topology::VirtualCore* getCollectorVirtualCore() const;
    const std::vector<mammut::topology::VirtualCore*>& getWorkersVirtualCore() const;

    const std::vector<mammut::topology::VirtualCore*>& getActiveVirtualCores() const;
    const std::vector<mammut::topology::VirtualCore*>& getUnusedVirtualCores() const;
private:
    /**
     * Generates mapping indexes. They are indexes to be used on
     * _availableVirtualCores vector to get the corresponding virtual core
     * where a specific node must be mapped.
     * @param emitterIndex The index of the emitter.
     * @param firstWorkerIndex The index of the first worker (the others follow).
     * @param collectorIndex The index of the collector (if present).
     */
    void getMappingIndexes(size_t& emitterIndex,
                           size_t& firstWorkerIndex,
                           size_t& collectorIndex);

    /**
     * Computes the mapping order of virtual cores for linear
     * mapping.
     */
    void computeVcOrderLinear();

    /**
     * Performs a linear mapping of the nodes on the available virtual cores.
     */
    void performLinearMapping();

    KnobConfMapping _confMapping;
    KnobConfSNodeMapping _confEmitterMapping;
    KnobConfSNodeMapping _confCollectorMapping;
    KnobConfHyperthreading _confHyperthreading;

    // The available virtual cores, sorted according to the mapping strategy.
    std::vector<mammut::topology::VirtualCore*>& _vcOrder;
    std::vector<mammut::topology::VirtualCore*> _vcOrderLinear;
    std::vector<mammut::topology::VirtualCore*> _vcOrderCacheEfficient;

    std::vector<mammut::topology::VirtualCore*> _activeVirtualCores;
    std::vector<mammut::topology::VirtualCore*> _unusedVirtualCores;

    mammut::topology::VirtualCore* _emitterVirtualCore;
    mammut::topology::VirtualCore* _collectorVirtualCore;
    std::vector<mammut::topology::VirtualCore*> _workersVirtualCores;
    AdaptiveNode* _emitter;
    AdaptiveNode* _collector;
    const KnobWorkers& _knobWorkers;
    mammut::topology::Topology* _topologyHandler;
};

class KnobFrequency: public Knob{
public:
    KnobFrequency(KnobConfFrequencies confFrequency,
                  const mammut::Mammut& mammut,
                  bool useTurboBoost,
                  StrategyUnusedVirtualCores unusedVc,
                  const KnobMapping& knobMapping);
    bool needsCalibration() const;
    void changeValueReal(double v);
    std::vector<double> getAllowedValues() const;
private:
    void applyUnusedVCStrategySame(const std::vector<mammut::topology::VirtualCore*>& unusedVc, mammut::cpufreq::Frequency v);
    void applyUnusedVCStrategyOff(const std::vector<mammut::topology::VirtualCore*>& unusedVc);
    void applyUnusedVCStrategyLowestFreq(const std::vector<mammut::topology::VirtualCore*>& unusedVc);
    void applyUnusedVCStrategy(mammut::cpufreq::Frequency v);

    KnobConfFrequencies _confFrequency;
    std::vector<double> _allowedValues;
    mammut::cpufreq::CpuFreq* _frequencyHandler;
    mammut::topology::Topology* _topologyHandler;
    mammut::cpufreq::CpuFreq* _cpufreqHandle;
    StrategyUnusedVirtualCores _unusedVc;
    const KnobMapping& _knobMapping;
};


/**
 * Knobs values can be:
 *  - Real: e.g. for workers it will get value between 1 and the maximum number of cores.
 *  - Relative: They will always assume values in the range [0.0, 100.0]
 */
typedef enum{
    KNOB_VALUE_UNDEF = 0,
    KNOB_VALUE_REAL,
    KNOB_VALUE_RELATIVE
}KnobValueType;

std::string knobTypeToString(KnobType kv);

class KnobsValues{
private:
    KnobValueType _type;
    double _values[KNOB_TYPE_NUM];
public:
    KnobsValues(KnobValueType type = KNOB_VALUE_UNDEF):_type(type){
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            _values[i] = 0;
        }
    }

    void swap(KnobsValues& x){
        using std::swap;

        swap(_type, x._type);
        swap(_values, x._values);
    }

    inline KnobsValues(const KnobsValues& other){
        _type = other._type;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            _values[i] = other._values[i];
        }
    }

    inline KnobsValues& operator=(KnobsValues other){
        swap(other);
        return *this;
    }

    inline bool areRelative() const{return _type == KNOB_VALUE_RELATIVE;}

    inline bool areReal() const{return _type == KNOB_VALUE_REAL;}

    inline double& operator[](KnobType idx){
        return _values[idx];
    }

    inline double operator[](KnobType idx) const{
        return _values[idx];
    }
};


inline std::ostream& operator<<(std::ostream& os, const KnobsValues& obj){
    os << "[";
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        os << obj[(KnobType) i] << ", ";
    }
    os << "]";
    return os;
}

inline bool operator==(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(lhs[(KnobType) i] !=
           rhs[(KnobType) i]){
            return false;
        }
    }
    return true;
}

inline bool operator!=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator==(lhs,rhs);
}

inline bool operator<(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
        if(lhs[(KnobType) i] <
           rhs[(KnobType) i]){
            return true;
        }else if(lhs[(KnobType) i] >
                 rhs[(KnobType) i]){
            return false;
        }
    }
    return false;
}

inline bool operator>(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    return operator< (rhs,lhs);
}

inline bool operator<=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator> (lhs,rhs);
}

inline bool operator>=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator< (lhs,rhs);
}

}

#endif /* SRC_KNOB_HPP_ */
