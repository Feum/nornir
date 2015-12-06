/*
 * configuration.hpp
 *
 * Created on: 05/12/2015
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

#ifndef ADAPTIVE_FASTFLOW_CONFIGURATION_HPP_
#define ADAPTIVE_FASTFLOW_CONFIGURATION_HPP_

#include "knob.hpp"
#include "node.hpp"
#include "parameters.hpp"

#include <mammut/utils.hpp>

namespace adpff{

class FarmConfiguration: public mammut::utils::NonCopyable {
private:
    Knob* _knobs[KNOB_TYPE_NUM];
    const Parameters& _p;
    std::vector<KnobsValues> _combinations;
    void combinations(std::vector<std::vector<double> > array, size_t i,
                      std::vector<double> accum);
    void setRelativeValues(const KnobsValues& values);
    void setRealValues(const KnobsValues& values);
public:
    template <typename lb_t, typename gt_t>
    FarmConfiguration(const Parameters& p, ff::ff_farm<lb_t, gt_t>& farm):_p(p){
        AdaptiveNode* emitter = dynamic_cast<AdaptiveNode*>(farm.getEmitter());
        AdaptiveNode* collector = dynamic_cast<AdaptiveNode*>(farm.getCollector());
        std::vector<AdaptiveNode*> w = convertWorkers(farm.getWorkers());
        ff::ff_gatherer* gt = farm.getgt();

        _knobs[KNOB_TYPE_WORKERS] = new KnobWorkers(p.knobWorkers, emitter,
                                                    collector, gt, w);
        _knobs[KNOB_TYPE_MAPPING] = new KnobMapping(p.knobMapping,
                                                    p.knobMappingEmitter,
                                                    p.knobMappingCollector,
                                                    p.knobHyperthreading,
                                                    p.mammut,
                                                    emitter,
                                                    collector,
                                                    *((KnobWorkers*)_knobs[KNOB_TYPE_WORKERS]));
        _knobs[KNOB_TYPE_FREQUENCY] = new KnobFrequency(p.knobFrequencies,
                                                        p.mammut,
                                                        p.turboBoost,
                                                        p.strategyUnusedVirtualCores,
                                                        *((KnobMapping*)_knobs[KNOB_TYPE_MAPPING]));

        std::vector<std::vector<double>> values;
        std::vector<double> accum;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            values.push_back(_knobs[i]->getAllowedValues());
        }
        combinations(values, 0, accum);
    }


    ~FarmConfiguration();

    /**
     * Gets all the possible combinations of knobs values.
     * @return A std::vector containing all the possible combinations
     *         of knobs values.
     */
    const std::vector<KnobsValues>& getAllRealCombinations() const;

    /**
     * Sets the highest frequency to reduce the reconfiguration time.
     */
    void setFastReconfiguration();

    /**
     * Returns a specified knob.
     * @param t The type of the knob to return.
     * @return The specified knob.
     */
    const Knob* getKnob(KnobType t) const;

    /**
     * Sets all the knobs to their maximum.
     */
    void maxAllKnobs();

    /**
     * Returns the real value of a specific knob.
     * @param t The type of the knob.
     * @return The real value of the specified knob.
     */
    double getRealValue(KnobType t) const;

    /**
     * Returns the real values for all the knobs.
     * @return The real values for all the knobs.
     */
    KnobsValues getRealValues() const;

    /**
     * Sets values for the knobs (may be relative or real).
     * @param values The values of the knobs.
     */
    void setValues(const KnobsValues& values);

};

}
#endif /* ADAPTIVE_FASTFLOW_CONFIGURATION_HPP_ */
