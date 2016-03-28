/*
 * mdfi.cpp
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

#include "mdfi.hpp"

namespace nornir{
namespace dataflow{

std::vector<OutputToken>* Mdfi::compute(){
    StreamElem* args[MAX_INSTRUCTION_INPUTS];
    for(uint i=0; i<dInput; i++)
        args[i]=tInput[i].getTask();
    StreamElem** result=comp->compute(args);
    std::vector<OutputToken> *v = new std::vector<OutputToken>(dOutput);
    for(uint i=0; i<dOutput; i++)
        /*Set output tokens.*/
        (*v)[i]=OutputToken(result[i],dest[i]);
    if(result!=NULL && result!=args)
        delete[] result;
    return v;
}

void Mdfi::updateDestinations(int* v){
    for(uint i=0; i<dOutput; i++){
        if(!dest[i].isOutStream()&&!dest[i].isNull())
            dest[i].setMdfiId(v[dest[i].getMdfId()]);
    }
}

void Mdfi::reset(unsigned long int newId){
    setGid(newId);
    for(uint i=0; i<dInput; i++)
        tInput[i].clear();
}

}
}


