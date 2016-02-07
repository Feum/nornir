/*
 * manager.hpp
 *
 * Created on: 23/03/2015
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

/*! \file manager.hpp
 * \brief Implementation of an adaptive fastflow farm.
 *
 * To let an existing fastflow farm-based adaptive, follow these steps:
 *  1. Emitter, Workers and Collector of the farm must extend
 *     adpff::AdaptiveNode instead of ff::ff_node
 *  2. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyWorkersChange virtual
 *     method.
 */

#ifndef ADAPTIVE_FASTFLOW_FARM_HPP_
#define ADAPTIVE_FASTFLOW_FARM_HPP_

//#if !defined(TRACE_FASTFLOW) || !defined(FF_TASK_CALLBACK)
//#error "TRACE_FASTFLOW and FF_TASK_CALLBACK macros need to be defined."
//#endif

#include "configuration.hpp"
#include "knob.hpp"
#include "parameters.hpp"
#include "predictors.hpp"
#include "node.hpp"
#include "utils.hpp"

#include <ff/farm.hpp>
#include <mammut/module.hpp>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>

#include <cmath>
#include <iostream>
#include <limits>

namespace adpff{

class Parameters;

//TODO REMOVE USING
using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

struct MonitoredSample;

/*!
 * \class ManagerFarm
 * \brief This class manages the adaptivity in farm based computations.
 *
 * This class manages the adaptivity in farm based computations.
 */
template <typename lb_t = ff::ff_loadbalancer, typename gt_t = ff::ff_gatherer>
class ManagerFarm: public Thread{
    friend class PredictorSimple;
    friend class PredictorLinearRegression;
    friend class RegressionData;
    friend class RegressionDataServiceTime;
    friend class RegressionDataPower;
    friend class Calibrator;
    friend class CalibratorLowDiscrepancy;
public:
    /**
     * Creates a farm adaptivity manager.
     * @param farm The farm to be managed.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFarm(ff_farm<lb_t, gt_t>* farm, Parameters adaptivityParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFarm();

    /**
     * Function executed by this thread.
     */
    void run();
private:
    // The managed farm.
    ff_farm<lb_t, gt_t>* _farm;

    // Flag for checking farm termination.
    volatile bool _terminated;

    // The parameters used to take management decisions.
    Parameters _p;

    // Starting time of the manager.
    double _startTimeMs;

    // The cpufreq module.
    CpuFreq* _cpufreq;

    // The energy counter.
    Counter* _counter;

    // The task module.
    TasksManager* _task;

    // The topology module.
    Topology* _topology;

    // The emitter (if present).
    AdaptiveNode* _emitter;

    // The collector (if present).
    AdaptiveNode* _collector;

    // The vector of active workers.
    std::vector<AdaptiveNode*> _activeWorkers;

    // Monitored samples;
    Smoother<MonitoredSample>* _samples;

    // The current configuration of the farm.
    FarmConfiguration _configuration;

    // The number of tasks processed since the last reconfiguration.
    double _totalTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the number of tasks
    // that still needs to be processed by the application.
    uint64_t _remainingTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the deadline of
    // the application.
    time_t _deadline;

    // Milliseconds timestamp of the last store of a sample.
    double _lastStoredSampleMs;

    // The calibrator of the predictors.
    Calibrator* _calibrator;

#ifdef DEBUG_MANAGER
    ofstream samplesFile;
#endif

    /**
     * Set a specified domain to the highest frequency.
     * @param domain The domain.
     */
    void setDomainToHighestFrequency(const Domain* domain);

    /**
     * Returns the primary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The primary value of a sample according to
     * the required contract.
     */
    double getPrimaryValue(const MonitoredSample& sample) const;

    /**
     * Returns the secondary value of a sample according to
     * the required contract.
     * @param sample The sample.
     * @return The secondary value of a sample according to
     * the required contract.
     */
    double getSecondaryValue(const MonitoredSample& sample) const;
    /**
     * Returns the primary value according to the required contract.
     * @return The primary value according to the required contract.
     */
    double getPrimaryValue() const;

    /**
     * Returns the secondary value according to the required contract.
     * @return The secondary value according to the required contract.
     */
    double getSecondaryValue() const;

    /**
     * Changes the knobs.
     */
    void changeKnobs();

    /**
     * Send data to observer.
     **/
    void observe();

    /**
     * Asks the workers for their samples.
     */
    void askForWorkersSamples();

    /**
     * Obtain workers samples.
     * @param sample A worker sample. It will be filled by this call with the
     *               global data of the farm.
     */
    void getWorkersSamples(WorkerSample& sample);

    /**
     * Resets the workers samples.
     */
    void resetSample();

    /**
     * Gets the consumed joules since the last reset and 
     * resets the counter.
     * @return The joules consumed since the last reset.
     */
    Joules getAndResetJoules();

    /**
     * Store a new sample.
     **/
    void storeNewSample();

    /**
     * Updates the tasks count.
     * @param sample The workers sample to be used for the update.
     */
    void updateTasksCount(WorkerSample& sample);

    /**
     * Returns true if the manager doesn't have still to check for a new
     * configuration.
     * @return True if the manager doesn't have still to check for a new
     * configuration.
     */
    bool persist() const;

    /**
     * Initializes the predictors
     */
    void initPredictors();

    /**
     * Operations that need to take place before running the nodes.
     */
    void initNodesPreRun();

    /**
     * Operations that need to take place after running the nodes.
     */
    void initNodesPostRun();

    /**
     * Cleans the adaptive nodes.
     */
    void cleanNodes();

    /**
     * Initializes the samples.
     * return A samples smoother with no recorded samples.
     */
    Smoother<MonitoredSample>* initSamples() const;
};

}

#include "manager.cpp"

#endif /* ADAPTIVE_FASTFLOW_FARM_HPP_ */
