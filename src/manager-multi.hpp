/*
 * manager-multi.hpp
 *
 * Created on: 22/07/2016
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

#ifndef NORNIR_MANAGER_MULTI_HPP_
#define NORNIR_MANAGER_MULTI_HPP_

#include "manager.hpp"

namespace nornir{

class ManagerMulti: public mammut::utils::Thread{
private:
    mammut::Mammut _m;
    mammut::topology::Topology* _topology;
    mammut::cpufreq::CpuFreq* _cpufreq;
    ff::SWSR_Ptr_Buffer _qIn;
    ff::SWSR_Ptr_Buffer _qOut;
    std::vector<Manager*> _activeManagers;
    std::vector<mammut::topology::PhysicalCoreId> _allCores;
    std::map<Manager*, std::vector<mammut::topology::PhysicalCoreId> > _allocatedCores;
    std::map<Manager*, std::vector<KnobsValues> > _allocations;
    std::vector<std::vector<size_t> > _allocationsCombinations;

    void inhibitAll(Manager* except);

    void disinhibitAll();

    std::vector<mammut::topology::PhysicalCoreId> getAvailablePhysicalCores() const;

    std::map<KnobsValues, double> invertMap(const std::map<KnobsValues, double>& map) const;

    void updateAllocations(Manager* m);

    void updateAllocations();

    void combinations(std::vector<std::vector<size_t> > array, size_t i, std::vector<size_t> accum);

    /**
     * Finds the best allocation.
     * @return A vector of indexes (one per manager). Each index is the position of
     * the chosen allocation inside the corresponding allocations vector.
     */
    std::vector<size_t> findBestAllocation();

    void applyNewAllocation();

    void applyWattsCorrection();
public:
    ManagerMulti();

    void addManager(Manager* m);

    /**
     * Returns a terminated manager if present, NULL otherwise.
     * @return A terminated manager if present, NULL otherwise.
     */
    Manager* getTerminatedManager();

    void run();
};

}

#endif /* NORNIR_MANAGER_MULTI_HPP_ */
