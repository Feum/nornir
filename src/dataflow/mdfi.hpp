/*
 * mdfi.hpp
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

#ifndef NORNIR_DF_MDFI_HPP_
#define NORNIR_DF_MDFI_HPP_

#include "skeleton/computable.hpp"
#include "tokens.hpp"

#include <cassert>
#include <iostream>
#include <vector>

namespace nornir{
namespace dataflow{

#define MAX_INSTRUCTION_INPUTS 100

/**Macro data flow instruction.**/
class Mdfi{
private:
    InputToken *tInput; ///<Input tokens. Are numbered from 0 to \e dInput -1.
    TokenId *dest; ///<Destinations of the output. Are numbered from 0 to \e dOutput -1.
    Computable* comp; ///<\e Computable to compute.
    unsigned int dInput, ///<Input size.
        dOutput, ///<Output size.
        id; ///<Instruction's identifier.
    unsigned long int graphId; ///<Id of the graph wich belongs the instruction.
public:
    /**
     * Constructor of the instruction.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param i Instruction's identifier.
     * \param nInput Size of the input.
     * \param nOutput Size of the output.
     */
    inline Mdfi(Computable* c, unsigned int i, unsigned int nInput,unsigned int nOutput):
    comp(c),dInput(nInput),dOutput(nOutput),id(i),graphId(0){
        tInput=new InputToken[nInput];
        dest=new TokenId[nOutput];
        assert(dInput <= MAX_INSTRUCTION_INPUTS);
    }

    /**
     * Copy constructor.
     * \param ins Reference to instruction to copy.
     */
    inline Mdfi(const Mdfi& ins):
            dest(ins.dest),comp(ins.comp),dInput(ins.dInput),dOutput(ins.dOutput),id(ins.id),graphId(ins.graphId){
        tInput=new InputToken[dInput];
        dest=new TokenId[dOutput];
        for(uint i=0; i<dOutput; i++)
            dest[i]=ins.dest[i];
        assert(dInput <= MAX_INSTRUCTION_INPUTS);
    }

    /**
     * Destructor of the instruction.
     */
    inline ~Mdfi(){
        delete[] tInput;
        delete[] dest;
    }

    /**
     * Computes the result of the instruction.
     * \return A vector of output tokens.
     */
    std::vector<OutputToken>* compute();

    /**
     * Checks if the instruction is fireable.
     * \return \e true if the instruction is fireable, otherwise returns \e false.
     */
    inline bool isFireable(){
        uint i=0;
        while(i<dInput){
            if(!tInput[i].isPresent()) return false;
            i++;
        }
        return true;
    }

    /**
     * Sets an input task.
     * \param t The task to set.
     * \param i The position where to set the task.
     * \return \e false if \e i is higher than \e dInput, otherwise returns false.
     */
    inline bool setInput(StreamElem* t,uint i){
        if(i<dInput) {
            tInput[i].setTask(t);
            return true;
        }else
            return false;
    }

    /**
     * Sets the \e i-th destination. Overwrites it if is present.
     * \param i The destination index.
     * \param d The destination to set.
     */
    inline void setDestination(unsigned int i,const TokenId& d){
        dest[i]=d;
    }

    /**
     * Returns the id of the instruction.
     * \return The id of the instruction.
     */
    inline unsigned int getId(){
        return id;
    }

    /**
     * Sets the id of the instruction.
     * \param i The id of the instruction.
     */
    inline void setId(unsigned int i){
        id=i;
    }

    /**
     * Returns the size of the input.
     * \return The size of the input.
     */
    inline unsigned int getInputSize(){
        return dInput;
    }

    /**
     * Returns the size of the output.
     * \return The size of the output.
     */
    inline unsigned int getOutputSize(){
        return dOutput;
    }

    /**
     * Returns the id of the graph.
     * \return The id of the graph.
     */
    inline unsigned long int getGid(){
        return graphId;
    }

    /**
     * Sets the id of the graph.
     * \param newd The new id of the graph.
     */
    inline void setGid(unsigned long int newd){
        graphId=newd;
        for(uint i=0; i<dOutput; i++)
            dest[i].setGraphId(newd);
    }

    /**
     * Returns a pointer to the \e Computable that the instruction have to compute.
     * \return A pointer to the \e Computable that the instruction have to compute.
     */
    inline Computable* getComputable(){
        return comp;
    }

    /**
     * Updates the destinations of the instructions.
     * \param v \e v[i] is the new identifier of the instruction that previously has \e i as identifier.
     */
    void updateDestinations(int* v);

    /**
     * Resets the instruction.
     * \param newId The new id of the graph wich belongs the instruction.
     */
    void reset(unsigned long int newId);
};

}
}

#endif /* NORNIR_DF_MDFI_HPP_ */

