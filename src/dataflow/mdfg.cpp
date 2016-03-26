/*
 * mdfg.cpp
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

#include "mdfg.hpp"

namespace nornir{
namespace dataflow{

Mdfg::Mdfg():higherId(10),nextId(0),id(0){
    instructions=(Mdfi**)malloc(sizeof(Mdfi*)*(higherId));
}

Mdfg::Mdfg(Computable* c):higherId(1),nextId(1),id(0){
    instructions=(Mdfi**)malloc(sizeof(Mdfi*)*(higherId));
    instructions[0]=new Mdfi(c,0,1,1);
    TokenId d;
    d.setOutputStream();
    /**Set the output stream as instruction's output.**/
    instructions[0]->setDestination(0,d);
}

Mdfg::Mdfg(Computable* c, int dInput, int dOutput):higherId(1),nextId(1),id(0){
    instructions=(Mdfi**)malloc(sizeof(Mdfi*)*(higherId));
    instructions[0]=new Mdfi(c,0,dInput,dOutput);
}

Mdfg::Mdfg(const Mdfg& g,unsigned long int gid):higherId(g.higherId),nextId(g.nextId),id(gid){
    instructions=(Mdfi**)malloc(sizeof(Mdfi*)*(nextId));
    unsigned int i=0;
    while(i<nextId){
        instructions[i]=new Mdfi(*(g.instructions[i]));
        instructions[i]->setGid(gid);
        i++;
    }
}

Mdfg::~Mdfg(){
    unsigned int i=0;
    while(i<nextId){
        delete instructions[i];
        i++;
    }
    free(instructions);
}

void Mdfg::reset(unsigned long int newId){
    id=newId;
    for(unsigned int i=0; i<nextId; i++)
        instructions[i]->reset(newId);
}

}
}

