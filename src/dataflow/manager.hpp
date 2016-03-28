/*
 * manager.hpp
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

#ifndef NORNIR_DF_MANAGER_HPP
#define NORNIR_DF_MANAGER_HPP

#include "interpreter.hpp"
#include "tokens.hpp"
#include "mdfg.hpp"
#include "stream.hpp"
#include "hashMap.hpp"
#include "skeleton/computable.hpp"
#include "skeleton/farm.hpp"
#include "skeleton/pipeline.hpp"
#include "skeleton/ewc.hpp"
#include <typeinfo>
#include <errno.h>
#include <deque>
#include <vector>

namespace nornir{
namespace dataflow{

/**
 * This is a recursive function that compiles a skeleton into a macro data flow graph.
 * \param c The skeleton to compile.
 */
Mdfg* compile(Computable* c);

/**
 * The manager of muskel.
 * The manager executes the skeleton (or the macro data flow graph).
 */
class Manager{
private:
    InputStream *in;
    OutputStream *out;
    /**A pointer to the interpreter.**/
    Interpreter *intr;
    /**The macro data flow graph.**/
    Mdfg *graph;
    /**True if the graph is a compiled graph.**/
    const bool compiled;
    /**Next free graph identifier.**/
    unsigned long int nextGraphId;
#ifdef POOL
    /**Pool of usable graphs.**/
    std::deque<Mdfg*> *pool;
#ifndef MAXPOOLSIZE
#define MAXPOOLSIZE 25000
#endif
#endif
    /**Instances of the graph.**/
    hashMap<Mdfg*> *graphs;
    /**
     * Computed results. It's necessary to save them into a map for preserving
     * the order of the task received from the input stream. The results will
     * be periodically send to the output stream.
     **/
    hashMap<StreamElem*> *result;
    /**Fireable instructions.**/
    std::deque<Mdfi*> *fireable;

    /**
     * Number of instructions to send to the interpreter before waiting the
     * results.
     **/
    unsigned long int groupSize,
    /**Number of instructions sent to interpreter.**/
        executed,
    /**Number of results not yet calculated.**/
        taskSent,
    /**Index of the last task sent to the output stream.**/
        lastSent;
    StreamElem** tempTask;
#ifdef COMPUTE_COM_TIME
    unsigned long acc;
#endif
public:
    /**
     * Constructor of the manager.
     * \param i The input stream.
     * \param o The output stream.
     * \param parDegree Number of workers to activate.
     * \param groupSize Number of instructions to send to the interpreter
     *                  before waiting the results.
     * \param c The skeleton to be executed.
     */
    Manager(InputStream *i,OutputStream *o, int parDegree, unsigned long int groupSize, Computable *c);

    /**
     * Constructor of the manager.
     * \param i The input stream.
     * \param o The output stream.
     * \param parDegree Number of workers to activate.
     * \param groupSize Number of instructions to send to the interpreter
     *                  before waiting the results.
     * \param graph The macro data flow graph to be executed.
     */
    Manager(InputStream *i, OutputStream *o, int parDegree, unsigned long int groupSize, Mdfg *graph );

    /**
     * Denstructor of the manager.
     */
    ~Manager();

    /**
     * Prints some stats.
     * \param out The stream where print the stats.
     */
    void stats(std::ostream& out);

    /**
     * Executes the skeleton (or the macro data flow graph).
     */
    void exec();
private:
    /**
     * Tries to receive new tasks from the input stream.
     **/
    void getFromInput();

    /**
     * Sends new results to the output stream, preserving the order of the
     * task received from the input stream.
     */
    void flushOnStream();

};
}
}
#endif // NORNIR_DF_MANAGER_HPP
