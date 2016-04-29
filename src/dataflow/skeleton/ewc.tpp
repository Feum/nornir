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
void ReduceEmitter<T>::compute(void){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) receiveData();
    int dim=task->getSize();
    int mod=dim%chunkNum;
    int size=dim/chunkNum;
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
        sendData(toAdd, (Computable*) i);
    }
    if(autoDelete) delete task;
#endif
}


template <typename T, T*(*fun)(T*,T*)>
void ReduceWorker<T, fun>::compute(void){
    void* t = receiveData();
    ArrayWrapper<T*>* w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai=(ArrayIndexes<T*>*)t;
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    void* toDelete = t;
    w1=(ArrayWrapper<T*>*)t;
    nElem=w1->getSize();
    first=0;
#endif
    if(nElem-first==1){
#ifndef NOCOPY
        t = w1->get(first);
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
    t = x;
    delete toDelete;
#endif
    }
    sendData(t);
}


template<typename T, T*(*fun)(T*,T*)>
ReduceCollector<T, fun>::ReduceCollector(int cn):chunkNum(cn){;}

template<typename T, T*(*fun)(T*,T*)>
void ReduceCollector<T, fun>::compute(void){
    T *p,*q;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai;
    ai=(ArrayIndexes<T*>*)receiveData(0);
    p=ai->getArray()->get(ai->geti());
    delete ai;
    ai=(ArrayIndexes<T*>*)receiveData(1);
    q=ai->getArray()->get(ai->geti());
    delete ai;
#else
    p=(T*)receiveData(0);
    q=(T*)receiveData((Computable*) 1);
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
        a=(T*)receiveData((Computable*) i);
#endif
        x=fun(x,a);
    }
    sendData(x);
}

template <typename T>
MapEmitter<T>::MapEmitter(uint numWorkers, bool autoDelete):numWorkers(numWorkers),autoDelete(autoDelete){;}

template <typename T>
void MapEmitter<T>::compute(void){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) receiveData();
    uint dim=task->getSize();
    uint mod=dim%numWorkers;
    uint size=dim/numWorkers;
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
    ArrayWrapper<void*>* toAdd;
    for(uint i=0; i<numWorkers; i++){
        if(mod && i==numWorkers-mod){
            size+=1;
            mod=0;
        }
        toAdd=new ArrayWrapper<void*>(size);
        for(uint j=0; j<size; j++){
            toAdd->set(j,task->get(k));
            k++;
        }
        sendData(toAdd, (Computable*) i);
    }
    if(autoDelete) delete task;
#endif
}

template <typename T, typename V, V*(*fun)(T*) >
void MapWorker<T, V, fun>::compute(void){
    ArrayWrapper<void*> *w1;
    int nElem,first;
    void* t = receiveData();
#ifdef NOCOPY
    ArrayIndexes<void*>* ai=(ArrayIndexes<void*>*) t;
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    w1=(ArrayWrapper<void*>*) t;
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
    sendData(t);
}

template <typename V>
MapCollector<V>::MapCollector(uint nWorkers):nWorkers(nWorkers){;}

template <typename V>
void MapCollector<V>::compute(void){
#ifdef NOCOPY
    ArrayIndexes<void*>* ai;
    for(uint i=0; i<nWorkers-1; i++){
        ai=(ArrayIndexes<void*>*)receiveData((Computable*) i);
        delete ai;
    }
    ai = (ArrayIndexes<void*>*)receiveData((Computable*) (nWorkers-1));
    sendData(ai->getArray());
    delete ai;
#else
    uint size=0;
    for(uint i=0; i<nWorkers; i++)
        size+=((ArrayWrapper<void*>*) receiveData((Computable*) i))->getSize();


    ArrayWrapper<V*> *aw=new ArrayWrapper<V*>(size),*tempAw;
    uint tempSize,k=0;
    for(uint i=0; i<nWorkers; i++){
        tempAw=((ArrayWrapper<V*>*) receiveData((Computable*) i));
        tempSize=tempAw->getSize();
        for(uint j=0; j<tempSize; j++){
            aw->set(k,tempAw->get(j));
            ++k;
        }
        delete tempAw;
    }
    sendData(aw);
#endif
}

}
}
