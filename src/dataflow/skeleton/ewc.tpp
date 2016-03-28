/*
 * ewc.tpp
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

template <typename T>
ReduceEmitter<T>::ReduceEmitter(int cn, bool autoDelete):
    chunkNum(cn),autoDelete(autoDelete){;}

template <typename T>
StreamElem** ReduceEmitter<T>::compute(StreamElem** t){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) t[0];
    int dim=task->getSize();
    int mod=dim%chunkNum;
    int size=dim/chunkNum;
    StreamElem** toRet=new StreamElem*[chunkNum];
#ifdef NOCOPY
    int l,h=0;
    for(uint i=0; i<chunkNum; i++){
        l=h;
        if(mod && i==chunkNum-mod){
            size+=1;
            mod=0;
        }
        h=l+size;
        toRet[i]=new ArrayIndexes<T*>(task,l,h);
    }
#else
    ArrayWrapper<T*>* toAdd;
    int k=0;
    for(uint i=0; i<chunkNum; i++){
        if(mod && i==chunkNum-mod){
            size+=1;
            mod=0;
        }
        toAdd=new ArrayWrapper<T*>(size);
        for(int j=0; j<size; j++){
            toAdd->set(j,task->get(k));
            k++;
        }
        toRet[i]=toAdd;
    }
    if(autoDelete) delete task;
#endif
    return toRet;
}


template <typename T, T*(*fun)(T*,T*)>
StreamElem** ReduceWorker<T, fun>::compute(StreamElem** t){
    ArrayWrapper<T*>* w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai=(ArrayIndexes<T*>*)t[0];
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    StreamElem* toDelete=t[0];
    w1=(ArrayWrapper<T*>*)t[0];
    nElem=w1->getSize();
    first=0;
#endif
    if(nElem-first==1){
#ifndef NOCOPY
        t[0]=w1->get(first);
        delete toDelete;
#endif
    }else{
        T *a=w1->get(first),*b=w1->get(first+1);
        T *x=fun(a,b);
        for(int i=first+2; i<nElem;i++){
            a=w1->get(i);
            x=fun(x,a);
        }
#ifdef NOCOPY
    w1->set(first,x);
#else
    t[0]=x;
    delete toDelete;
#endif
    }
    return t;
}


template<typename T, T*(*fun)(T*,T*)>
ReduceCollector<T, fun>::ReduceCollector(int cn):chunkNum(cn){;}

template<typename T, T*(*fun)(T*,T*)>
StreamElem** ReduceCollector<T, fun>::compute(StreamElem** t){
    StreamElem** toRet=new StreamElem*[1];
    T *p,*q;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai;
    ai=(ArrayIndexes<T*>*)t[0];
    p=ai->getArray()->get(ai->geti());
    delete ai;
    ai=(ArrayIndexes<T*>*)t[1];
    q=ai->getArray()->get(ai->geti());
    delete ai;
#else
    p=(T*)t[0];
    q=(T*)t[1];
#endif
    T* x=fun(p,q);
    T* a;
    for(int i=2; i<chunkNum; i++){
#ifdef NOCOPY
        ArrayIndexes<T*>* ai;
        ai=(ArrayIndexes<T*>*)t[i];
        a=ai->getArray()->get(ai->geti());
        if(i==chunkNum-1) delete ai->getArray();
        delete ai;
#else
        a=(T*)t[i];
#endif
        x=fun(x,a);
    }
    toRet[0]=x;
    return toRet;
}

template <typename T>
MapEmitter<T>::MapEmitter(uint numWorkers, bool autoDelete):numWorkers(numWorkers),autoDelete(autoDelete){;}

template <typename T>
StreamElem** MapEmitter<T>::compute(StreamElem** t){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) t[0];
    uint dim=task->getSize();
    uint mod=dim%numWorkers;
    uint size=dim/numWorkers;
    StreamElem** toRet=new StreamElem*[numWorkers];
#ifdef NOCOPY
    int l,h=0;
    for(uint i=0; i<numWorkers; i++){
        l=h;
        if(mod && i==numWorkers-mod){
            size+=1;
            mod=0;
        }
        h=l+size;
        toRet[i]=new ArrayIndexes<T*>(task,l,h);
    }
#else
    int k=0;
    ArrayWrapper<StreamElem*>* toAdd;
    for(uint i=0; i<numWorkers; i++){
        if(mod && i==numWorkers-mod){
            size+=1;
            mod=0;
        }
        toAdd=new ArrayWrapper<StreamElem*>(size);
        for(uint j=0; j<size; j++){
            toAdd->set(j,task->get(k));
            k++;
        }
        toRet[i]=toAdd;
    }
    if(autoDelete) delete task;
#endif
    return toRet;
}

template <typename T, typename V, V*(*fun)(T*) >
StreamElem** MapWorker<T, V, fun>::compute(StreamElem** t){
    ArrayWrapper<StreamElem*> *w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<StreamElem*>* ai=(ArrayIndexes<StreamElem*>*)t[0];
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    w1=(ArrayWrapper<StreamElem*>*)t[0];
    nElem=w1->getSize();
    first=0;
#endif
    V* y;
    T* x;
    for(int i=first; i<nElem;i++){
        x=(T*)w1->get(i);
        y=fun(x);
        w1->set(i,y);
    }
    return t;
}

template <typename V>
MapCollector<V>::MapCollector(uint nWorkers):nWorkers(nWorkers){;}

template <typename V>
StreamElem** MapCollector<V>::compute(StreamElem** t){
    StreamElem** toRet=new StreamElem*[1];
#ifdef NOCOPY
    ArrayIndexes<StreamElem*>* ai;
    for(uint i=0; i<nWorkers-1; i++){
        ai=(ArrayIndexes<StreamElem*>*)t[i];
        delete ai;
    }
    ai=(ArrayIndexes<StreamElem*>*)t[nWorkers-1];
    toRet[0]=ai->getArray();
    delete ai;
    return toRet;
#else
    uint size=0;
    for(uint i=0; i<nWorkers; i++)
        size+=((ArrayWrapper<StreamElem*>*) t[i])->getSize();


    ArrayWrapper<V*> *aw=new ArrayWrapper<V*>(size),*tempAw;
    uint tempSize,k=0;
    for(uint i=0; i<nWorkers; i++){
        tempAw=((ArrayWrapper<V*>*) t[i]);
        tempSize=tempAw->getSize();
        for(uint j=0; j<tempSize; j++){
            aw->set(k,tempAw->get(j));
            ++k;
        }
        delete tempAw;
    }
    toRet[0]=aw;
    return toRet;
#endif
}

}
}
