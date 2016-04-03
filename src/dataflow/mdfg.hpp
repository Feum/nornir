/*
 * mdfg.hpp
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

#ifndef NORNIR_DF_MDFG_HPP_
#define NORNIR_DF_MDFG_HPP_

#include "mdfi.hpp"

namespace nornir{
namespace dataflow{

/**
 * Macro data flow graph.
 * Accessory methods are provided to interact with and use the graph.
 *
 * The Input Mdfi must have an unique input token.
 *
 * The Output Mdfi must have an unique destination.
 *
 * Usage example:
 * \code
 *      using namespace nornir;
 *      using namespace nornir::dataflow;
 *      Mdfg* graph = new Mdfg();
 *      uint i1, i2, i3, i4, i5, i6;
 *
 *      i1 = graph->createFirstMdfi(seq1, 2); <--Creates the first instruction that compute the sequential skeleton seq1
 *                                            The instruction have 2 output. (The number of input is implicit one)
 *      i2 = graph->createMdfi(seq2, 1, 1);
 *      i3 = graph->createMdfi(seq3, 1, 1);
 *      i4 = graph->createMdfi(seq4, 2, 1);
 *      i5 = graph->createMdfi(seq5, 1, 1);
 *      i6 = graph->createLastMdfi(seq6, 1);  <--Creates the last instruction. The instruction have 1 input. (The number
 *                                            of output is implicit one)
 * \endcode
 * To link them together:
 * \code
 *      graph->link(i1, 0, i2, 0);
 *      graph->link(i1, 1, i3, 0);
 *      graph->link(i2, 0, i4, 0);
 *      graph->link(i3, 0, i4, 1);
 *      graph->link(i4, 0, i5, 0);
 *      graph->link(i5, 0, i6, 0);//TODO AGGIUNGERE IMG GRAFO.
 * \endcode
 */
class Mdfg{
private:
    /**Next free index**/
    uint _nextId;
    /**Graph's id.**/
    ulong _id;
    /**Instructions of the graph.**/
    std::vector<Mdfi> _instructions;
public:
    /**
     * Constructor of the graph.
     */
    Mdfg();

    /**
     * Builds a graph with only one instruction.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     */
    Mdfg(Computable* c);

    /**
     * Builds a graph with only one instruction. It's used for the construction of the emitter and
     * the collector of a map/reduce.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     */
    Mdfg(Computable* c, int dInput, int dOutput);

    /**
     * Copy constructor.
     * \param g Reference to graph to copy
     * \param gid The new id of the graph's copy.
     */
    Mdfg(const Mdfg& g, ulong gid);

    /**
     * Graph's destructor.
     */
    ~Mdfg();


    inline ulong getId(){
        return _id;
    }


    /**
     * Returns the number of instructions.
     * \return Number of instructions.
     */
    inline uint getNumMdfi(){
        return _instructions.size();
    }

    /**
     * Creates the first instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param nOutput Size of the output.
     * \return The id of the created instruction.
     */
    inline uint createFirstMdfi(Computable* c, unsigned int nOutput){
        return createMdfi(c, 1, nOutput);
    }

    /**
     * Creates the first instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return The id of the created instruction.
     */
    inline uint createFirstMdfi(Mdfi* in){
        if(in->getInputSize() != 1){
            throw std::runtime_error("First Mdfi must have only 1 input.");
        }else{
            return createMdfi(in);
        }
    }

    /**
     * Creates an instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param nInput Size of the input.
     * \param nOutput Size of the output.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Computable* c, uint nInput, uint nOutput){
        _instructions.emplace_back(c, _nextId, nInput, nOutput);
        ++_nextId;
        return _instructions.size() - 1;
    }

    /**
     * Creates an instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Mdfi* in){
        _instructions.emplace_back(*in);
        _instructions.back().setId(_nextId);
        ++_nextId;
        return _instructions.size() - 1;
    }

    /**
     * Creates an instruction of the graph starting from an
     * instruction of another graph.
     * \param g The other graph.
     * \param id The id of the mdfi to be copied.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Mdfg* g, uint id){
        _instructions.emplace_back(g->_instructions.at(id));
        _instructions.back().setId(_nextId);
        ++_nextId;
        return _instructions.size() - 1;
    }

    /**
     * Creates the last instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param nInput Size of the input.
     * \return The id of the created instruction.
     */
    inline uint createLastMdfi(Computable* c, unsigned int nInput){
        uint torid = createMdfi(c, nInput, 1);
        TokenId d;
        d.setOutputStream();
        /**Set the output stream as instruction's output.**/
        _instructions.at(torid).setDestination(0, d);
        return torid;
    }

    /**
     * Creates the last instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return The id of the created instruction.
     */
    inline uint createLastMdfi(Mdfi* in){
        if(in->getOutputSize() != 1){
            throw std::runtime_error("Last Mdfi must have only 1 output.");
        }else{
            uint torid = createMdfi(in);
            TokenId d;
            d.setOutputStream();
            /**Set the output stream as instruction's output.**/
            _instructions.at(torid).setDestination(0, d);
            return torid;
        }
    }

    /**
     * Creates the last instruction of the graph starting
     * from the instruction of another graph.
     * \param g The other graph.
     * \param id The id of the instruction in the other graph.
     * \return The id of the created instruction.
     */
    inline uint createLastMdfi(Mdfg* g, uint id){
        if(g->_instructions.at(id).getOutputSize() != 1){
            throw std::runtime_error("Last Mdfi must have only 1 output.");
        }else{
            uint torid = createMdfi(g, id);
            TokenId d;
            d.setOutputStream();
            /**Set the output stream as instruction's output.**/
            _instructions.at(torid).setDestination(0, d);
            return torid;
        }
    }

    /**
     * Returns the id of the first instruction of the graph.
     * \return The id of the created instruction.
     */
    inline uint getFirst(){
        if(_instructions.size()){
            return 0;
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
    }

    /**
     * Returns the id of the last instruction of the graph.
     * \return The id of the created instruction.
     */
    inline uint getLast(){
        if(_instructions.size()){
            return _instructions.size() - 1;
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
    }

    /**
     * Links two instructions. The \e y-th output of the instruction \e x will be linked to
     * the \e w-th output of the instruction \e z.
     * e.g.:
     * \code
     * graph.link(2, 0, 3, 0); <--Link the first output of mdfi with id 2 to the first input of mdfi with id 3.
     * \endcode
     */
    inline bool link(uint x, uint y, uint z, uint w){
        TokenId tid(z, w);
        /**Check the indexes.**/
        if((y >= _instructions.at(x).getOutputSize()) ||
           (w >= _instructions.at(z).getInputSize())){
            return false;
        }else{
            _instructions.at(x).setDestination(y, tid);
            return true;
        }
    }

    /**
     * Updates the destinations of an instruction.
     * \param id The id of the Mdfi.
     * \param v The new destinations.
     */
    inline void updateDestinations(uint id, int* v){
        if(id >= _instructions.size()){
            throw std::runtime_error("Non existing mdfi.");
        }
        _instructions.at(id).updateDestinations(v);
    }

    inline bool isFireable(uint id) const{
        return _instructions.at(id).isFireable();
    }

    /**
     * ATTENTION: SHOULD ONLY BE CALLED AFTER THE GRAPH HAS BEEN CREATED.
     *            If instructions are added/removed, this pointer will be
     *            invalidated.
     */
    inline Mdfi* getMdfi(uint id){
        return &(_instructions.at(id));
    }

    /**
     * Resets all the instructions.
     * \param newId The new id of the graph.
     */
    void reset(ulong newId);

};
}
}

#endif /* MDFG_HPP_ */
