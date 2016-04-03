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
 *      Mdfi *i1,*i2,*i3,*i4,*i5,*i6;
 *
 *      i1=graph->createFirstMdfi(seq1,2); <--Creates the first instruction that compute the sequential skeleton seq1
 *                                            The instruction have 2 output. (The number of input is implicit one)
 *      i2=graph->createMdfi(seq2,1,1);
 *      i3=graph->createMdfi(seq3,1,1);
 *      i4=graph->createMdfi(seq4,2,1);
 *      i5=graph->createMdfi(seq5,1,1);
 *      i6=graph->createLastMdfi(seq6,1);  <--Creates the last instruction. The instruction have 1 input. (The number
 *                                            of output is implicit one)
 * \endcode
 * To link them together:
 * \code
 *      graph->link(i1,0,i2,0);
 *      graph->link(i1,1,i3,0);
 *      graph->link(i2,0,i4,0);
 *      graph->link(i3,0,i4,1);
 *      graph->link(i4,0,i5,0);
 *      graph->link(i5,0,i6,0);//TODO AGGIUNGERE IMG GRAFO.
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
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createFirstMdfi(Computable* c, unsigned int nOutput){
        return createMdfi(c, 1, nOutput);
    }

    /**
     * Creates the first instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createFirstMdfi(Mdfi* in){
        if(in->getInputSize() != 1){
            return NULL;
        }else{
            return createMdfi(in);
        }
    }

    /**
     * Creates an instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param nInput Size of the input.
     * \param nOutput Size of the output.
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createMdfi(Computable* c, uint nInput, uint nOutput){
        _instructions.emplace_back(c, _nextId, nInput, nOutput);
        ++_nextId;
        return &(_instructions.back());
    }

    /**
     * Creates an instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createMdfi(Mdfi* in){
        _instructions.emplace_back(*in);
        _instructions.back().setId(_nextId);
        ++_nextId;
        return &(_instructions.back());
    }

    /**
     * Creates the last instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param nInput Size of the input.
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createLastMdfi(Computable* c, unsigned int nInput){
        Mdfi* tor = createMdfi(c, nInput, 1);
        TokenId d;
        d.setOutputStream();
        /**Set the output stream as instruction's output.**/
        tor->setDestination(0, d);
        return tor;
    }

    /**
     * Creates the last instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return A pointer to the created instruction.
     */
    inline Mdfi* createLastMdfi(Mdfi* in){
        if(in->getOutputSize() != 1){
            return NULL;
        }else{
            Mdfi* tor = createMdfi(in);
            TokenId d;
            d.setOutputStream();
            /**Set the output stream as instruction's output.**/
            tor->setDestination(0, d);
            return tor;
        }
    }

    /**
     * Returns a pointer to the first instruction of the graph.
     * \return A pointer to the first instruction of the graph.
     */
    inline Mdfi* getFirst(){
        if(_instructions.size()){
            return &(_instructions.front());
        }else{
            return NULL;
        }
    }

    /**
     * Returns a pointer to the instruction with identifier equal to \e i.
     * \param i The identifier of the instruction.
     */
    inline Mdfi* getMdfi(uint i){
        return &(_instructions.at(i));
    }

    /**
     * Returns a pointer to the last instruction of the graph.
     * \return A pointer to the last instruction of the graph.
     */
    inline Mdfi* getLast(){
        if(_instructions.size()){
            return &(_instructions.back());
        }else{
            return NULL;
        }
    }

    /**
     * Links two instructions. The \e y-th output of the instruction \e x will be linked to
     * the \e w-th output of the instruction \e z.
     * e.g.:
     * \code
     * graph.link(x,0,z,0); <--Link the first output of x to the first input of z.
     * \endcode
     */
    inline bool link(Mdfi* x, uint y, Mdfi* z, uint w){
        TokenId tid(z->getId(), w);
        /**Check the indexes.**/
        if((y >= x->getOutputSize()) || (w >= z->getInputSize())){
            return false;
        }else{
            x->setDestination(y, tid);
            return true;
        }
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
