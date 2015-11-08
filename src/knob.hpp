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

class Knob: public mammut::utils::NonCopyable{
public:
    Knob():_relativeValue(-1), _realValue(-1){;}

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
     * Returns the current relative value [0, 100] of this knob.
     * @return The current relative value [0, 100] of this knob.
     */
    double getRelativeValue() const;

    /**
     * Returns the current real value of this knob.
     * @return The current real value of this knob.
     */
    double getRealValue() const;

    /**
     * Returns true if the best value for this knob needs to be automatically
     * found.
     * @return True if the best value for this knob needs to be automatically
     * found.
     */
    bool autoFind() const;

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

    double _relativeValue;
    double _realValue;
};

class KnobWorkers: public Knob{
public:
    KnobWorkers(KnobConfWorkers confWorkers, ff::ff_farm<>& farm);
    void changeValueReal(double v);
    std::vector<double> getAllowedValues() const;
    uint getNumActiveWorkers() const;
    const std::vector<AdaptiveNode*>& getActiveWorkers() const;
private:
    /**
     * Returns all the workers of the farm.
     * @param farm The farm.
     * @return All the workers of the farm
     */
    std::vector<AdaptiveNode*> getAllWorkers(const ff::ff_farm<>& farm) const;

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
    ff::ff_farm<>& _farm;
    AdaptiveNode* _emitter;
    AdaptiveNode* _collector;
    const std::vector<AdaptiveNode*> _allWorkers;
    std::vector<AdaptiveNode*> _activeWorkers;
    std::vector<double> _knobValues;
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
    void changeValueReal(double v);
    std::vector<double> getAllowedValues() const;
private:
    void applyUnusedVCStrategyOff(const std::vector<mammut::topology::VirtualCore*>& unusedVc);
    void applyUnusedVCStrategyLowestFreq(const std::vector<mammut::topology::VirtualCore*>& vc);
    void applyUnusedVCStrategy();

    KnobConfFrequencies _confFrequency;
    std::vector<double> _allowedValues;
    mammut::cpufreq::CpuFreq* _frequencyHandler;
    mammut::topology::Topology* _topologyHandler;
    mammut::cpufreq::CpuFreq* _cpufreqHandle;
    StrategyUnusedVirtualCores _unusedVc;
    const KnobMapping& _knobMapping;
};

}

#endif /* SRC_KNOB_HPP_ */
