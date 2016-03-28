/*
 * interpreter.cpp
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

#include "interpreter.hpp"

namespace nornir{
namespace dataflow{

WorkerMdf::WorkerMdf(ff::dynqueue* buffer):_buffer(buffer){;}

void WorkerMdf::compute(Mdfi* t){
    std::vector<OutputToken>* temp = t->compute();
    _buffer->push(temp);
}

Interpreter::Interpreter(size_t maxWorkers):_maxWorkers(maxWorkers){
    _p = new nornir::Parameters("parameters.xml", "archdata.xml");
    _o = new nornir::Observer;
    _p->observer = _o;

    _accelerator = new nornir::FarmAccelerator<Mdfi>(_p);

    /**Creates the SPSC queues.**/
    _buffers = new ff::dynqueue*[_maxWorkers];
    for(size_t i = 0; i < _maxWorkers; i++){
        _buffers[i] = new ff::dynqueue();
    }
    /**Adds the workers to the farm.**/
    for(size_t i = 0; i < _maxWorkers; ++i){
        _accelerator->addWorker(new WorkerMdf(_buffers[i]));
    }

    _accelerator->start();
}

Interpreter::~Interpreter(){
     _accelerator->wait();
     for(int i = 0; i < _maxWorkers; i++){
         delete _buffers[i];
     }
     delete[] _buffers;
     delete _p;
     delete _o;
}

/**
 * Waits for the results of the instructions passed to the interpreter.
 * \param r A queue of output tokens. In this queue the interpreter will puts the computed results.
 */
int Interpreter::wait(ff::squeue<OutputToken>& r){
    int acc = 0;
    std::vector<OutputToken>* temp;
    void* task;
    for(size_t i = 0; i < _maxWorkers; i++){
#ifdef COMPUTE_COM_TIME
        unsigned long t1;
        while(true){
            t1 = ff::getusec();
            if(!buffers[i]->pop(&task)) break;
            acc += ff::getusec()-t1;
#else
        while(_buffers[i]->pop(&task)){
#endif
            temp = (std::vector<OutputToken>*) task;
            for(int j = 0; j< (int) temp->size(); j++)
                r.push_back((*temp)[j]);
            delete temp;
        }
    }
    return acc;
}

}
}

