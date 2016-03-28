/*
 * computable.hpp
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

#ifndef NORNIR_DF_COMPUTABLE_HPP_
#define NORNIR_DF_COMPUTABLE_HPP_

#include "../stream.hpp"

namespace nornir{
namespace dataflow{

/**
 * \class Computable
 * This class represents a generic computation.
 * If a class extends \e Computable, its "compute" function can be computed sequentially in a Macro Data Flow Instruction.
 * \e Compute is a function that takes an array of Task* and returns an array of Task*. The array must be dinamically
 * allocated using new[].
 *
 * This is an example of a function that takes \e n values as input and returns \e m values.
 *
 * \code
 * #include "nornir/dataflow/skeleton/computable.hpp"
 *
 * class myCompute: public skel::Computable{
 * public:
 * Task** compute(Task** input){
 *      IntTask* x;
 *      //Creation of the array to return.
 *      //(To note that the array is created using new[].
 *      Task** toReturn=new Task*[n];
 *      for(int i=0; i<n; i++){
 *          //Gets the input task (to note the cast).
 *          x=(IntTask*) input[i];
 *          //Insertion of the elements to return.
 *          //toReturn[i] will be sent to the MDFi linked
 *          //to the i-th output.
 *          toReturn[i]=new IntArrayTask(x[i]->get()+10);
 *      }
 *      return toReturn;
 *  }
 *};
 * \endcode
 */

class Computable{
public:
    /**
     * Destructor of the computable.
     */
    virtual ~Computable(){;}

    /**
     * \param t Array of Task*.
     * \return Array of Task* as result.
     *
     * This method computes the result.
     *
     */
    virtual StreamElem** compute(StreamElem** t) = 0;
};

}
}

#endif /* NORNIR_DF_COMPUTABLE_HPP_ */
