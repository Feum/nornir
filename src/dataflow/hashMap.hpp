/*
 * hashMap.hpp
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

#ifndef NORNIR_DF_HASHMAP_HPP_
#define NORNIR_DF_HASHMAP_HPP_
#include <list>

namespace nornir{
namespace dataflow{

template <class T> class hashMap{
private:
    std::list<std::pair<long int,T> >** l;
public:
    hashMap();

    ~hashMap();

    void put(long int index, const T& x);

    bool get(long int index, T* p) const;

    void erase(long int index);

    bool getAndErase(long int index, T* p);
};

}
}

#include "hashMap.tpp"

#endif /**NORNIR_DFHASHMAP_HPP_**/
