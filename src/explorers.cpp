/*
 * explorers.cpp
 *
 * Created on: 10/07/2015
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

#include "explorers.hpp"
#include "configuration.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_EXPLORERS
#define DEBUG(x) do { std::cerr << "[Explorers] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{
    Explorer::Explorer(const FarmConfiguration& configuration):_configuration(configuration){;}

    ExplorerRandom::ExplorerRandom(const FarmConfiguration& configuration):
            Explorer(configuration){
        srand(time(NULL));
    }

    KnobsValues ExplorerRandom::nextRelativeKnobsValues() const{
        KnobsValues r(KNOB_VALUE_RELATIVE);
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            if(_configuration.getKnob((KnobType) i)->needsCalibration()){
                r[(KnobType)i] = rand() % 100;
            }else{
                /**
                 * If we do not need to automatically find the value for this knob,
                 * then it has only 0 or 1 possible value. Accordingly, we can set
                 * it to any value. Here we set it to 100.0 for readability.
                 */
                r[(KnobType)i] = 100.0;
            }
        }
        return r;
    }

    void ExplorerRandom::reset(){;}

    ExplorerLowDiscrepancy::ExplorerLowDiscrepancy(
                const FarmConfiguration& configuration,
                StrategyExploration explorationStrategy):
            Explorer(configuration), _explorationStrategy(explorationStrategy){
        uint d = 0;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            if(_configuration.getKnob((KnobType) i)->needsCalibration()){
                ++d;
            }
        }
        DEBUG("[Explorator] Will generate low discrepancy points in " << d <<
              " dimensions");
        const gsl_qrng_type* generatorType;
        switch(_explorationStrategy){
            case STRATEGY_EXPLORATION_NIEDERREITER:{
                generatorType = gsl_qrng_niederreiter_2;
            }break;
            case STRATEGY_EXPLORATION_SOBOL:{
                generatorType = gsl_qrng_sobol;
            }break;
            case STRATEGY_EXPLORATION_HALTON:{
                generatorType = gsl_qrng_halton;
            }break;
            case STRATEGY_EXPLORATION_HALTON_REVERSE:{
                generatorType = gsl_qrng_reversehalton;
            }break;
            default:{
                throw std::runtime_error("ExplorerLowDiscrepancy: Unknown "
                                         "exploration strategy: " +
                                         _explorationStrategy);
            }break;
        }
        _generator = gsl_qrng_alloc(generatorType, d);
        _normalizedPoint = new double[d];
    }

    ExplorerLowDiscrepancy::~ExplorerLowDiscrepancy(){
        gsl_qrng_free(_generator);
        delete[] _normalizedPoint;
    }


    KnobsValues ExplorerLowDiscrepancy::nextRelativeKnobsValues() const{
        KnobsValues r(KNOB_VALUE_RELATIVE);
        gsl_qrng_get(_generator, _normalizedPoint);
        size_t nextCoordinate = 0;
        for(size_t i = 0; i < KNOB_TYPE_NUM; i++){
            if(_configuration.getKnob((KnobType) i)->needsCalibration()){
                r[(KnobType)i] = _normalizedPoint[nextCoordinate]*100.0;
                ++nextCoordinate;
            }else{
                /**
                 * If we do not need to automatically find the value for this knob,
                 * then it has only 0 or 1 possible value. Accordingly, we can set
                 * it to any value. Here we set it to 100.0 for readability.
                 */
                r[(KnobType)i] = 100.0;
            }
        }
        return r;
    }

    void ExplorerLowDiscrepancy::reset(){
        gsl_qrng_init(_generator);
    }
}
