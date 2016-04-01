/*
 * hashMap.tpp
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

namespace nornir{
namespace dataflow{

#ifndef HASHSIZE
#define HASHSIZE 100000
#endif

template <class T>
hashMap<T>::hashMap(){
    l=new std::list<std::pair<long int,T> >*[HASHSIZE];
    for(int i=0; i<HASHSIZE; i++)
        l[i]=new std::list<std::pair<long int,T> >();
}

template <class T>
hashMap<T>::~hashMap(){
    for(int i=0; i<HASHSIZE; i++)
        delete l[i];
    delete[] l;
}

template <class T>
void hashMap<T>::put(long int index, const T& x){
    l[index%HASHSIZE]->push_back(std::make_pair(index,x));
}

template <class T>
bool hashMap<T>::get(long int index, T* p) const{
    std::list<std::pair<long int,T> >* v = l[index%HASHSIZE];
    typename std::list<std::pair<long int,T> >::iterator it=v->begin();
    while(it!=v->end() && (*it).first!=index)
        ++it;

    if(it!=v->end()){
        *p=(*it).second;
        return true;
    }
    else return false;
}

template <class T>
void hashMap<T>::erase(long int index){
    std::list<std::pair<long int,T> >* v=l[index%HASHSIZE];
    typename std::list<std::pair<long int,T> >::iterator it=v->begin();
    while(it!=v->end() && (*it).first!=index)
        ++it;
    if(it!=v->end()) v->erase(it);
}

template <class T>
bool hashMap<T>::getAndErase(long int index, T* p){
    std::list<std::pair<long int,T> >* v=l[index%HASHSIZE];
    typename std::list<std::pair<long int,T> >::iterator it=v->begin();
    while(it!=v->end() && (*it).first!=index)
        ++it;

    if(it!=v->end()){
        *p=(*it).second;
        v->erase(it);
        return true;
    }
    else return false;
}

}
}
