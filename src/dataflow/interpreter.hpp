/*
 * interpreter.hpp
 *
 * Created on: 26/03/2016
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

#ifndef NORNIR_DF_INTERPRETER_HPP_
#define NORNIR_DF_INTERPRETER_HPP_

#include "mdfi.hpp"
#include "../interface.hpp"
#include "../external/fastflow/ff/squeue.hpp"
#include <vector>
#include <queue>

namespace nornir{
namespace dataflow{

/**
 * This is the worker of the fastflow's farm.
 */
class WorkerMdf: public nornir::Worker<Mdfi>{
private:
    ff::dynqueue* _buffer;
public:
    WorkerMdf(ff::dynqueue* buffer);

    /**
     * Computes a macro data flow instruction.
     * \param task A macro data flow instruction.
     * \return The result of the computation.
     */
    void compute(Mdfi* t);
};

/**
 * This is the interpreter of the macro data flow instructions.
 */
class Interpreter{
private:
    ff::dynqueue** _buffers;
    nornir::FarmAccelerator<Mdfi>* accelerator;
    int parDegree;
    nornir::Parameters* p;
    nornir::Observer* o;
public:
    /**
     * Constructor of the interpreter.
     * \param parDegree Parallelism degree (number of worker to activate).
     */
    Interpreter(int parDegree);

    /**
     * Destructor of the interpreter.
     */
    ~Interpreter();

    /**
     * Computes asynchronously a macro data flow instruction.
     * \param instr The instruction to compute.
     */
    inline void exec(Mdfi* instr){
        accelerator->offload(instr);
    }

    /**
     * Prints some stats.
     * \param out The stream where print the stats.
     */
    inline void stats(std::ostream& out=std::cout){
        out << "Stats not available." << std::endl;
    }

    /**
     * Waits for the results of the instructions passed to the interpreter.
     * \param r A queue of output tokens. In this queue the interpreter will puts the computed results.
     */
    int wait(ff::squeue<OutputToken>& r);

    /**Stops the farm.**/
    inline void stop(){
        accelerator->offload(NULL);
    }
};

}
}

#endif /* NORNIR_DF_INTERPRETER_HPP_ */
