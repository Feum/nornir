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

EmitterWorkerCollector::EmitterWorkerCollector(Computable* emitter,
        Computable* worker, Computable* collector, uint nWorkers,
        bool deleteAll):
                nWorkers(nWorkers), emitter(emitter),worker(worker),
                collector(collector),deleteAll(deleteAll){;}

/**
 * Destructor of the EmitterWorkerCollector.
 * If deleteAll is true, deletes emitter, worker and collector.
 */
EmitterWorkerCollector::~EmitterWorkerCollector(){
    if(deleteAll){
        delete emitter;
        delete worker;
        delete collector;
    }
}

Task** EmitterWorkerCollector::compute(Task** t){
    Task **emitterResult,**fromWorker,
        **result=new Task*[nWorkers],*toWorker[1];
    emitterResult=emitter->compute(t);

    for(uint i=0; i<nWorkers; i++){
        toWorker[0]=emitterResult[i];
        fromWorker=worker->compute(toWorker);
        result[i]=fromWorker[0];
    }
    delete[] emitterResult;
    Task** fromCollector=collector->compute(result);
    delete[] result;
    return fromCollector;
}


template <typename T>
ReduceEmitter<T>::ReduceEmitter(int cn, bool autoDelete):
    chunkNum(cn),autoDelete(autoDelete){;}

template <typename T>
Task** ReduceEmitter<T>::compute(Task** t){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) t[0];
    int dim=task->getSize();
    int mod=dim%chunkNum;
    int size=dim/chunkNum;
    Task** toRet=new Task*[chunkNum];
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
Task** ReduceWorker<T, fun>::compute(Task** t){
    ArrayWrapper<T*>* w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai=(ArrayIndexes<T*>*)t[0];
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    Task* toDelete=t[0];
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
Task** ReduceCollector<T, fun>::compute(Task** t){
    Task** toRet=new Task*[1];
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
Task** MapEmitter<T>::compute(Task** t){
    ArrayWrapper<T*>* task=(ArrayWrapper<T*>*) t[0];
    uint dim=task->getSize();
    uint mod=dim%numWorkers;
    uint size=dim/numWorkers;
    Task** toRet=new Task*[numWorkers];
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
    ArrayWrapper<Task*>* toAdd;
    for(uint i=0; i<numWorkers; i++){
        if(mod && i==numWorkers-mod){
            size+=1;
            mod=0;
        }
        toAdd=new ArrayWrapper<Task*>(size);
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
Task** MapWorker<T, V, fun>::compute(Task** t){
    ArrayWrapper<Task*> *w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<Task*>* ai=(ArrayIndexes<Task*>*)t[0];
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    w1=(ArrayWrapper<Task*>*)t[0];
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
Task** MapCollector<V>::compute(Task** t){
    Task** toRet=new Task*[1];
#ifdef NOCOPY
    ArrayIndexes<Task*>* ai;
    for(uint i=0; i<nWorkers-1; i++){
        ai=(ArrayIndexes<Task*>*)t[i];
        delete ai;
    }
    ai=(ArrayIndexes<Task*>*)t[nWorkers-1];
    toRet[0]=ai->getArray();
    delete ai;
    return toRet;
#else
    uint size=0;
    for(uint i=0; i<nWorkers; i++)
        size+=((ArrayWrapper<Task*>*) t[i])->getSize();


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
