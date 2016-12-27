/*
 * manager.hpp
 *
 * Created on: 23/03/2015
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
 * @file manager.hpp
 * @brief Implementation of an adaptive fastflow farm.
 */

#ifndef NORNIR_FARM_HPP_
#define NORNIR_FARM_HPP_

#include "configuration.hpp"
#include "ffincs.hpp"
#include "knob.hpp"
#include "parameters.hpp"
#include "selectors.hpp"
#include "node.hpp"
#include "utils.hpp"

#include "external/mammut/mammut/module.hpp"
#include "external/mammut/mammut/utils.hpp"
#include "external/mammut/mammut/mammut.hpp"
#include "external/orlog/src/orlog.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace nornir{

class Parameters;
class ManagerMulti;

//TODO REMOVE USING
using namespace std;
using namespace ff;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

struct MonitoredSample;

typedef struct{
    size_t numSteps;
    double bandwidthAccuracy;
    double powerAccuracy;
    double currentBandwidth;
    double currentPower;
    KnobsValues foundConfiguration;
    std::vector<double> performanceErrors;
    std::vector<double> powerErrors;
}SimulationResult;

// How to react to the calibration of other applications
// in order to do not interfere too much with them.
typedef enum{
    // Do nothing
    CALIBRATION_SHRINK_NONE = 0, //TODO: New managers arrived to the multimanager may not have space to calibrate.
    // Move all the threads on the same core.
    CALIBRATION_SHRINK_AGGREGATE,
    // Move all the threads on the same core and pause the process.
    CALIBRATION_SHRINK_PAUSE,
}CalibrationShrink;

/*!
 * \class Manager
 * \brief This class manages the adaptivity parallel applications.
 *
 * This class manages the adaptivity in parallel applications.
 */
class Manager: public Thread{
    friend class ManagerMulti;
public:
    Manager(Parameters adaptivityParameters);

    virtual ~Manager();
    /**
     * Function executed by this thread.
     * ATTENTION: The user must not call this one but 'start()'.
     */
    void run();

    /** 
     * Forces the manager to terminate.
     **/
    void terminate();
protected:
    // Flag for checking farm termination.
    volatile bool _terminated;

    // The parameters used to take management decisions.
    Parameters _p;

    // The cpufreq module.
    CpuFreq* _cpufreq;

    // The energy counter.
    Counter* _counter;

    // The task module.
    TasksManager* _task;

    // The topology module.
    Topology* _topology;

    // Monitored samples;
    Smoother<MonitoredSample>* _samples;

    // Variations
    Smoother<double>* _variations;

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

    // Inhibition flag.
    bool _inhibited;

    // The current configuration of the application.
    Configuration* _configuration;

    // The configuration selector.
    Selector* _selector;

    // Pid of the monitored process.
    pid_t _pid;
#ifdef DEBUG_MANAGER
    ofstream samplesFile;
#endif

    /**
     * Wait for the application to start and
     * sets the ProcessHandler to KnobMappingExternal if necessary.
     */
    virtual void waitForStart() = 0;

    /**
     * Asks the application for a sample.
     */
    virtual void askSample() = 0;

    /**
     * Obtain application sample.
     * @param sample An application sample. It will be filled by this call.
     */
    virtual void getSample(orlog::ApplicationSample& sample) = 0;

    /**
     * Manages a configuration change.
     */
    virtual void manageConfigurationChange() = 0;

    /**
     * Cleaning after termination.
     */
    virtual void clean() = 0;

    /**
     * Returns the execution time of the application (milliseconds).
     * - For ManagerFarm, it will be the execution time of the farm.
     * - For ManagerExternal, it will be the execution time between the construction
     * of the Application object and the call of the 'terminate' function on it
     * (application side).
     * - For ManagerBlackBox, it will be the execution time from the point where
     *   the manager is created to the end of the application.
     */
    virtual ulong getExecutionTime() = 0;

    /**
     * Updates the required bandwidth.
     */
    void updateRequiredBandwidth();

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
     * Returns true if the manager doesn't have still to check for a new
     * configuration.
     * @return True if the manager doesn't have still to check for a new
     * configuration.
     */
    bool persist() const;

    /**
     * Locks the knobs according to the selector/predictor.
     */
    void lockKnobs() const;

    /**
     * Initializes the samples.
     * return A samples smoother with no recorded samples.
     */
    Smoother<MonitoredSample>* initSamples() const;

    /**
     * Updates the tasks count.
     * @param sample The workers sample to be used for the update.
     */
    void updateTasksCount(orlog::ApplicationSample& sample);

    /**
     * Observes.
     **/
    void observe();

    /**
     * Decides the next configuration.
     */
    KnobsValues decide();

    /**
     * Move to the specified configuration.
     * @param kv The configuration
     */
    void act(KnobsValues kv, bool force = false);

    /**
     * Gets the consumed joules since the last reset and
     * resets the counter.
     * @return The joules consumed since the last reset.
     */
    Joules getAndResetJoules();

    /**
     * Logs the last observation.
     **/
    void logObservation();

    /**
     * Creates the selector.
     */
    Selector* createSelector() const;
private:
    void inhibit();

    void disinhibit();

    void shrink(CalibrationShrink type);

    void stretch(CalibrationShrink type);

    virtual void shrinkPause() = 0;

    virtual void stretchPause() = 0;

    /**
     * Updates the prediction models by letting them now
     * that there is an interference caused by another application.
     */
    void updateModelsInterference();


    /**
     * Wait until the models have been updated and they are
     * ready to be used.
     */
    void waitModelsInterferenceUpdate();

    /**
     * Returns the vector of physical cores used by the manager.
     * If the manager is still calibrating, it waits until
     * the manager finishes the calibration.
     *
     * @return The vector of physical cores used by the manager.
     */
    std::vector<PhysicalCoreId> getUsedCores();

    void allowCores(std::vector<mammut::topology::VirtualCoreId> ids);
};

class ManagerExternal: public Manager{
private:
    orlog::Monitor _monitor;
    void shrinkPause();
    void stretchPause();
public:
    /**
     * Creates an adaptivity manager for an external INSTRUMENTED application.
     * @param orlogChannel The name of the Orlog channel.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerExternal(const std::string& orlogChannel,
                    Parameters adaptivityParameters);

    /**
     * Creates an adaptivity manager for an external application.
     * @param orlogSocket The Orlog socket.
     * @param chid The channel id.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerExternal(nn::socket& orlogSocket,
                    int chid,
                    Parameters adaptivityParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerExternal();
protected:
    void waitForStart();
    void askSample();
    void getSample(orlog::ApplicationSample& sample);
    void manageConfigurationChange();
    void clean();
    ulong getExecutionTime();
};

class ManagerBlackBox: public Manager{
private:
    mammut::task::ProcessHandler* _process;
    double _startTime;
    double _lastTime;
    void shrinkPause();
    void stretchPause();
public:
    /**
     * Creates an adaptivity manager for an external NON-INSTRUMENTED application.
     * @param pid The identifier of an already running process.
     * @param adaptivityParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerBlackBox(pid_t pid, Parameters adaptivityParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerBlackBox();

    /**
     * Returns the pid of the monitored process.
     * @return The pid of the monitored process.
     */
    pid_t getPid() const;
protected:
    void waitForStart();
    void askSample();
    void getSample(orlog::ApplicationSample& sample);
    void manageConfigurationChange();
    void clean();
    ulong getExecutionTime();
};



/*!
 * \class ManagerFarm
 * \brief This class manages the adaptivity in farm based computations.
 *
 * This class manages the adaptivity in farm based computations.
 */
template <typename lb_t = ff::ff_loadbalancer, typename gt_t = ff::ff_gatherer>
class ManagerFarm: public Manager{
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
     * Simulates the execution.
     * ATTENTION: This is only meant to be used by developers.
     * @param configurationData The lines contained in configurationData file.
     * @param maxIterations The maximum number of iterations to be performed
     *        during calibration phase. If 0, there is no bound on the maximum
     *        number of iterations.
     */
    SimulationResult simulate(std::vector<std::string>& configurationData,
                              volatile bool* terminate,
                              size_t maxIterations = 0);
private:
    // The managed farm.
    ff_farm<lb_t, gt_t>* _farm;

    // The emitter (if present).
    AdaptiveNode* _emitter;

    // The collector (if present).
    AdaptiveNode* _collector;

    // The vector of active workers.
    std::vector<AdaptiveNode*> _activeWorkers;

    void waitForStart();
    void askSample();
    void getSample(orlog::ApplicationSample& sample);
    void initNodesPreRun();
    void initNodesPostRun();
    void manageConfigurationChange();
    void clean();
    ulong getExecutionTime();
    void shrinkPause();
    void stretchPause();
};

}

#include "manager.tpp"

#endif /* NORNIR_FARM_HPP_ */
