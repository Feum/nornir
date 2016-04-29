/*
 * tokens.hpp
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

#ifndef NORNIR_DF_TOKENS_HPP_
#define NORNIR_DF_TOKENS_HPP_

#include "stream.hpp"
#include <limits>

namespace nornir{
namespace dataflow{

#define MAXUNSLONGINT __LONG_MAX__ * 2UL + 1
#define MAXUNSINT __INT_MAX__ * 2U + 1

/**
 * A token's destination identifier is a triple <Graph identifier, Instruction identifier, Input slot identifier>
 */
class Mdfi;

class TokenId{
    friend class Mdfi;
private:
    unsigned long int graphId;///<The id of the graph.
    unsigned int mdfiId,///<The id of the instruction.
        tokenId;///<The id of the input slot.
    bool isOutputStream;///<True if the destination is the output stream.
public:
    /**
     * Constructor of the identifier.
     */
    inline TokenId():graphId(MAXUNSLONGINT),mdfiId(MAXUNSINT),
    tokenId(MAXUNSINT),isOutputStream(false){;}

    /**
     * Constructor of the identifier.
     * \param gid Identifier of the graph.
     * \param instrid Identifier of the instruction.
     * \param tokid Identifier of the input slot.
     */
    inline TokenId(unsigned long int gid, unsigned int instrid, unsigned int tokid):
            graphId(gid),mdfiId(instrid),tokenId(tokid),isOutputStream(false){;}

    /**
     * Constructor of the identifier.
     * \param instrid Identifier of the instruction.
     * \param tokid Identifier of the input slot.
     */
    inline TokenId(unsigned int instrid,unsigned int tokid):
            graphId(MAXUNSLONGINT),mdfiId(instrid),tokenId(tokid),isOutputStream(false){;}

    /**Sets the destination to the output stream.**/
    inline void setOutputStream(){
        isOutputStream=true;
    }

    /**Unsets the destination to the output stream.**/
    inline void unsetOutputStream(){
        isOutputStream=false;
    }

    /**
     * Checks if the destination is present.
     * \return \e True if the destination isn't present, \e false otherwise.
     */
    inline bool const isNull(){
        return graphId==MAXUNSLONGINT  && mdfiId==MAXUNSINT
                && tokenId==MAXUNSINT;
    }

    /**
     * Checks if the destination is the output stream.
     * \return \e True if the destination is the output stream, \e false otherwise.
     */
    inline bool isOutStream(){
        return isOutputStream;
    }

    /**
     * Returns the identifier of the destination graph.
     * \return The identifier of the destination graph.
     */
    inline unsigned long int getGraphId(){
        return graphId;
    }

    /**
     * Returns the identifier of the destination instruction.
     * \return The identifier of the destination instruction.
     */
    inline unsigned int getMdfId(){
        return mdfiId;
    }

    /**
     * Returns the destination instruction's slot.
     * \return The destination instruction's slot.
     */
    inline unsigned int getTokId(){
        return tokenId;
    }

    /**
     * Sets the identifier of the destination graph.
     * \param id The identifier of the destination graph.
     */
    inline void setGraphId(unsigned long int id){
        graphId=id;
    }

    /**
     * Sets the identifier of the destination instruction.
     * \param id The identifier of the destination instruction.
     */
    inline void setMdfiId(unsigned int id){
        mdfiId=id;
    }

    /**
     * Sets the identifier of the token.
     * \return The identifier of the token.
     */
    inline void setTokId(unsigned int id){
        tokenId=id;
    }
};

/**
 * An input token.
 */
class InputToken{
private:
    friend class Mdfi;
    void* task;///<The input task.
    bool presence;///<\e True if the task is present, \e false otherwise.
public:
    /**
     * Costructor of the token.
     * \param t The input task.
     */
    inline InputToken(void* t):task(t),presence(true){;}

    /**
     * Constructor of the token.
     */
    inline InputToken():task(0),presence(false){;}

    /**
     * Checks the presence of the task.
     * \return \e True if the task is present, otherwise returns \e false.
     */
    inline bool isPresent() const{
        return presence;
    }

    /**
     * Sets the input task.
     * \param t The task to set.
     */
    inline void setTask(void* t){
        task=t;
        presence=true;
    }

    /**
     * Returns the input task.
     * \return The input task.
     */
    inline void* getTask(){
        return task;
    }

    /**
     * Resets this input task.
     */
    inline void clear(){
        task=NULL;
        presence=false;
    }

};

/**
 * An output token is a pair <Task calculated,Destination of the token>.
 */
class OutputToken{
private:
    friend class Mdfi;
    void* result;
    TokenId dest;
public:
    /**Default constructor.**/
    inline OutputToken():result(NULL), dest(TokenId()){;}

    /**
     * Constructor of the token.
     * \param r A pointer to the computed task.
     * \param d The destination of the token.
     */
    inline OutputToken(void* r, const TokenId& d):result(r),dest(d){;}

    inline void clear(){result = NULL; dest = TokenId();}

    /**
     * Returns a pointer to the computed task.
     * \return A pointer to the computed task.
     */
    inline void* getResult(){return result;}

    /**
     * Returns the destination of the token.
     * \return The destination of the token.
     */
    inline TokenId getDest(){return dest;}
};

}
}
#endif /* NORNIR_DF_TOKENS_HPP_ */
