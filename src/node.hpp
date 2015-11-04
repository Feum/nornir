/*
 * node.hpp
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

#ifndef ADAPTIVE_FASTFLOW_NODE_HPP_
#define ADAPTIVE_FASTFLOW_NODE_HPP_

#include "./parameters.hpp"

#include <ff/node.hpp>
#include <mammut/mammut.hpp>

#undef DEBUG
#undef DEBUGB


namespace adpff{


/*!
 * \internal
 * \struct NodeSample
 * \brief Represents a sample of values taken from an adaptive node.
 *
 * This struct represents a sample of values taken from an adaptive node.
 */
typedef struct WorkerSample{
    // The percentage of time that the node spent on svc().
    double loadPercentage;

    // The number of computed tasks.
    double tasksCount;

    // The average service time (in ticks).
    double latency;

    // The bandwidth of the node.
    double bandwidthTotal;

    WorkerSample():loadPercentage(0), tasksCount(0),
                   latency(0), bandwidthTotal(0){;}

    WorkerSample& operator+=(const WorkerSample& rhs){
        loadPercentage += rhs.loadPercentage;
        tasksCount += rhs.tasksCount;
        latency += rhs.latency;
        bandwidthTotal += rhs.bandwidthTotal;
        return *this;
    }
}NodeSample;

inline WorkerSample operator+(WorkerSample lhs, const WorkerSample& rhs){
    lhs += rhs;
    return lhs;
}

/*!
 * \internal
 * \class ManagementRequestType
 * \brief Possible request types that a manager can make.
 */
typedef enum{
    // Get the current sample and reset it.
    MGMT_REQ_GET_AND_RESET_SAMPLE = 0,

    // Freezes the farm.
    MGMT_REQ_FREEZE,

    // Thaws the farm.
    MGMT_REQ_THAW
}ManagementRequestType;

/**
 * A management request.
 */
typedef struct{
    ManagementRequestType type;
    void* mark;
    size_t numWorkers;
}ManagementRequest;

/*!private
 * \class adpff_node
 * \brief This class wraps a ff_node to let it reconfigurable.
 *
 * This class wraps a ff_node to let it reconfigurable.
 */
class AdaptiveNode: public ff::ff_node{
private:
    friend class ManagerFarm;
    friend class KnobWorkers;
    friend class KnobMapping;

    volatile bool _started;
    volatile bool _terminated;
    volatile bool _goingToFreeze;
    mammut::task::TasksManager* _tasksManager;
    mammut::task::ThreadHandler* _thread;
    ManagementRequest _managementRequest;
    WorkerSample _sampleResponse;
    double _ticksPerNs;
    ticks _startTicks;
    ticks _ticksTot;
    ticks _tasksCount;

    // Queue used by the manager to notify that a request is present
    // on _managementRequest.
    ff::SWSR_Ptr_Buffer _managementQ;

    // Queue used by the node to notify that a response is present
    // on _sampleResponse.
    ff::SWSR_Ptr_Buffer _responseQ;

    /**
     * Converts ticks to seconds.
     * @param ticks The ticks to be converted.
     * @return The amount of seconds corresponding to the value of ticks.
     */
    double ticksToSeconds(double ticks) const;

    /**
     * Initializes the node.
     * @param mammut A Mammut handle.
     * @param ticksPerNs The number of ticks in a nanosecond.
     */
    void init(mammut::Mammut& mammut, double ticksPerNs);

    /**
     * Cleans then node.
     */
    void clean();

    /**
     * Moves this node on a specific virtual core.
     * @param vc The virtual core where this nodes must be moved.
     */
    void move(mammut::topology::VirtualCore* vc);

    /**
     * The result of askForSample call.
     * @param sample The statistics computed since the last
     *               time 'askForSample' has been called.
     * @param strategyPolling Strategy to apply if the queue is empty.
     * @param avgLatency Current average latency of the workers (in ns).
     */
    void getSampleResponse(WorkerSample& sample,
                           StrategyPolling strategyPolling,
                           double avgLatency);

    /**
    * Ask the node for a sample of the statistics computed since the last
    * time this method has been called.
    * The result can be retrieved with getSampleResponse call.
    */
    void askForSample();

    /**
     * Tells the node to freeze the farm.
     */
    void freezeAll(void* mark);

    /**
     * Thaws the farm.
     */
    void thawAll(size_t numWorkers);
    /**
     * Called on the node before starting it.
     * It must be called when the node is running.
     */
    void prepareToFreeze();

    /**
     * Called on the node before starting them.
     * It must be called when the node is frozen.
     */
    void prepareToRun();

    /**
     * Returns true if the farm has been terminated by the application,
     * false otherwise.
     * @return true if the farm has been terminated by the application,
     * false otherwise.
     */
    bool isTerminated() const;

    /**
     * Stores a sample.
     */
    void storeSample();

    /**
     * The callback that will be executed by the ff_node before
     * reading a task from the queue.
     * @param p Is the lb_t or gt_t in case of emitter or collector.
     */
    void callbackIn(void *p) CX11_KEYWORD(final);

    /**
     * Notify the end of stream.
     */
    void eosnotify(ssize_t) CX11_KEYWORD(final);

public:
    /**
     * Builds an adaptive node.
     */
    AdaptiveNode();

    /**
     * Destroyes this adaptive node.
     */
    ~AdaptiveNode();

    /**
     * This method can be implemented by the nodes to be aware of a change
     * in the number of workers.
     * When the farm is stopped and before running it again with the new
     * number of workers, this method is called. It is called on the emitter
     * (if present), on the collector (if present) and on all the workers of
     * the new configuration.
     * In this way, if needed action may be taken to prepare for the new
     * configuration (e.g. shared state modification, etc..).
     * @param oldNumWorkers The old number of workers.
     * @param newNumWorkers The new number of workers.
     */
    virtual void notifyWorkersChange(size_t oldNumWorkers,
                                     size_t newNumWorkers);
};

}

#endif /* ADAPTIVE_FASTFLOW_NODE_HPP_ */
