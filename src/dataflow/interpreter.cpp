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
#include <map>

namespace nornir{
namespace dataflow{

/**
 * This is a recursive function that compiles a skeleton into a macro data flow graph.
 * \param c The skeleton to compile.
 */
Mdfg* compile(Computable* c){
    const std::type_info& tc = typeid(*c);

    /**
     * Farm compiling. Creates a graph containing the instructions
     * of the worker.
     **/
    if(tc==typeid(Farm)){
        Farm* f=(Farm*) c;
        return compile(f->getWorker());
    /**
     * Pipeline compiling. Adds all istructions of the second stage in the
     * first stage of the pipeline and links it together.
     **/
    }else if(tc==typeid(Pipeline)){
        Pipeline* p = (Pipeline*) c;
        /**Compiles the two stages.**/
        Mdfg* g1 = compile(p->getFirstStage());
        Mdfg* g2 = compile(p->getSecondStage());
        bool lastSet = false;
        uint n = g2->getNumMdfi();
        uint added,
             oldSecondStageFirst, ///<The first instruction of the second stage of the pipeline.
             oldLast; ///<The last instruction of the first stage of the pipeline.

        /**
         * \e newIds is an array of Mdfi identifiers. \e newIds[i] is the new
         * identifier of the instruction of the second stage, that previously
         * has \e i as identifier.
         **/
        int *newIds = new int[n];
        oldLast = g1->getLast();
        oldSecondStageFirst = g2->getFirst();
        /**
         * Discriminates the case where the second stage has only one
         * instruction.
         **/
        if(n == 1){
            oldSecondStageFirst = g1->createLastMdfi(g2, 0);
            lastSet = true;
        }else{
            oldSecondStageFirst = g1->createMdfi(g2, 0);
        }
        newIds[0] = oldSecondStageFirst;
        /**
         * Adds the instructions of the second stage (except first and last
         * instruction).
         **/
        for(uint i = 1; i < n - 1; i++){
            added = g1->createMdfi(g2, i);
            newIds[i] = added;
        }

        /**If the last instruction isn't set.**/
        if(!lastSet){
            newIds[n - 1] = g1->createLastMdfi(g2, g2->getLast());
        }
        /**Updates destinations.**/
        for(uint i = 0; i < n; i++){
            g1->updateDestinations(newIds[i], newIds);
        }

        /**Links the two stages.**/
        g1->link(oldLast, 0, oldSecondStageFirst, 0);

        delete[] newIds;
        /**
         * Deletes the old graph of the second stage because all its
         * instructions were copied in the graph of the first stage.
         **/
        delete g2;
        return g1;
    /**Map and reduce compiling.**/
    }else if(tc==typeid(EmitterWorkerCollector)){
        EmitterWorkerCollector* ewc = (EmitterWorkerCollector*) c;
        int workersNum = ewc->getNWorkers();
        /**Compiles worker.**/
        Mdfg* worker = compile(ewc->getWorker());
        Mdfg* emitter = new Mdfg(ewc->getEmitter(), 1, workersNum);
        Mdfg* collector = new Mdfg(ewc->getCollector(), workersNum, 1);

        int workerSize = worker->getNumMdfi();
        /**
         * \e newWorkerIds is an array of Mdfi identifiers. \e newWorkerIds[i]
         * is the new identifier of the instruction of the worker, that
         * previously has \e i as identifier.
         **/
        int *newWorkerIds = new int[workerSize];
        /**
         * \e firstWorkerInstr and \e lastWorkerInstrs contain the first and
         * the last instructions of the workers.
         **/
        int *firstWorkerInstr = new int[workersNum];
        int *lastWorkerInstr = new int[workersNum];
        /**Adds the workers.**/
        for(int i = 0; i < workersNum; i++){
            for(int j = 0; j < workerSize; j++){
                newWorkerIds[j] = emitter->createMdfi(worker, j);
            }
            firstWorkerInstr[i] = newWorkerIds[0];
            lastWorkerInstr[i] = newWorkerIds[workerSize - 1];
            /**Updates the destinations.**/
            for(int j = 0; j < workerSize; j++){
                emitter->updateDestinations(newWorkerIds[j], newWorkerIds);
            }
        }
        delete[] newWorkerIds;
        /**Links the emitter to the workers.**/
        for(int i = 0; i < workersNum; i++){
            emitter->link(emitter->getFirst(), i, firstWorkerInstr[i], 0);
        }
        delete[] firstWorkerInstr;
        /**Adds the collector.**/
        uint firstCollInstr = emitter->createLastMdfi(collector, 0);

        /**Links the workers to collector.**/
        for(int i = 0; i < workersNum; i++){
            emitter->link(lastWorkerInstr[i], 0, firstCollInstr, i);
        }
        delete[] lastWorkerInstr;
        /**
         * Deletes the worker and the collector because all their instructions
         *  were copied in the graph of the emitter.
         **/
        delete worker;
        delete collector;

        return emitter;
    /**If c is an unknown skeleton, then it is a sequential skeleton.**/
    }else{
        return new Mdfg(c);
    }
}

class Scheduler: public nornir::Scheduler<Mdfi>{
private:
    InputStream *_in;
    OutputStream *_out;
    /**The macro data flow graph.**/
    Mdfg *_graph;
    /**True if the graph is a compiled graph.**/
    const bool _compiled;
    /**Next free graph identifier.**/
    unsigned long int _nextGraphId;
    /**Pool of usable graphs.**/
    std::deque<Mdfg*> *_pool;
    /**Instances of the graph.**/
    std::map<ulong, Mdfg*> *_graphs;
    /**
     * Computed results. It's necessary to save them into a map for preserving
     * the order of the task received from the input stream. The results will
     * be periodically send to the output stream.
     **/
    std::map<ulong, StreamElem*> *_result;

    unsigned long int
    /**Number of results not yet calculated.**/
        _taskSent,
    /**Index of the last task sent to the output stream.**/
        _lastSent;
    size_t _maxWorkers;
    size_t _numWorkers;
    size_t _groupSize;
    bool _orderedTasks;
#ifdef COMPUTE_COM_TIME
    unsigned long _acc;
#endif

    ff::uSWSR_Ptr_Buffer** _buffers;
    size_t _lastRcvId;
    ulong _tasksInside;
    ulong _graphsInside;

    inline void sendToWorkers(Mdfi* instr){
        ++_tasksInside;
        if(_orderedTasks){
            sendTo(instr, instr->getId() % _numWorkers);
        }else{
            send(instr);
        }
    }

    /**
     * - When we insert a graph we increase both graphsInside and tasksInside.
     * - When we pop a task the following situations may occur:
     *   - Decrease tasksInside and graphsInside.
     *   - Keep graphsInside untouched and change tasksInside (increase o
     *     decrease by n according to the specific graph structure).
     */

    inline bool keepInsertingGraphs() const{
        if(_groupSize){
            return _graphsInside < _groupSize;
        }else{
            return true;
        }
    }

    inline bool keepPoppingTasks() const{
        if(_groupSize){
            return _graphsInside > _groupSize;
        }else{
            return false;
        }
    }

    inline void updateCompleted(){
        OutputToken ot;
        ulong graphId;
        TokenId dest;
        Mdfi* ins;
        void* task;
        uint collected = 0;
        uint collectedTotal = 0;
        size_t startId = 0;
        size_t workerId = 0;
        startId = rand();
        do{
            collected  = 0;
            for(size_t i = 0; i < _maxWorkers; i++){
                workerId = (i + startId) % _maxWorkers;
#ifdef COMPUTE_COM_TIME
                unsigned long t1;
                if(true){
                    t1 = ff::getusec();
                    if(!buffers[workerId]->pop(&task)) break;
                    _acc += ff::getusec()-t1;
#else
                if(_buffers[workerId]->pop(&task)){
#endif
                    --_tasksInside;
                    ++collected;
                    ++collectedTotal;
                    for(uint j = 0; j < ((Mdfi*) task)->getNumOutTokens(); j++){
                        ot = *(((Mdfi*) task)->getOutToken(j));
                        dest = ot.getDest();
                        graphId = dest.getGraphId();
                        Mdfg* currentGraph = NULL;
                        /**
                         * If was the last instruction of a graph's copy, puts it
                         * into the output vector and delete
                         * the copy of the graph.
                         */
                        if(dest.isOutStream()){
                            assert(_result->emplace(std::piecewise_construct,
                                   std::forward_as_tuple(graphId),
                                   std::forward_as_tuple(ot.getResult())).second);

                            auto it = _graphs->find(graphId);
                            if(it != _graphs->end()){
                                currentGraph = it->second;
                                _graphs->erase(it);
                            }else{
                                throw std::runtime_error("Graph not found.");
                            }
#ifdef POOL
                            if(_pool->size() < MAXPOOLSIZE){
                                _pool->push_back(currentGraph);
                            }else{
                                delete currentGraph;
                            }
#else
                            delete currentGraph;
#endif
                            --_taskSent;
                            --_graphsInside;
                        }else{
                            /**Takes the pointer to the copy of the graph.**/
                            currentGraph = _graphs->at(graphId);
                            /**Takes the pointer to the instruction.**/
                            ins = currentGraph->getMdfi(dest.getMdfId());
                            /**Updates the instruction adding the input token.**/
                            ins->setInput(ot.getResult(), dest.getTokId());
                            /**
                             * If the instruction is fireable, adds it to the pool of
                             * fireable instructions.
                             **/
                            if(ins->isFireable()){
                                sendToWorkers(ins);
                            }
                        }
                    }
                }
            }
        }while(collected && keepPoppingTasks());
    }


    void getFromInput(StreamElem* next){
        Mdfg *newGraph;
        Mdfi *first;
#ifdef POOL
        if(_pool->size()){
            newGraph = _pool->front();
            _pool->pop_front();
            newGraph->reset(_nextGraphId);
            assert(_graphs->insert(std::pair<ulong, Mdfg*>(_nextGraphId, newGraph)).second);
        }else{
            newGraph = new Mdfg(*_graph, _nextGraphId);
            assert(_graphs->insert(std::pair<ulong, Mdfg*>(_nextGraphId, newGraph)).second);
        }
#else
        newGraph = new Mdfg(*_graph, _nextGraphId);
        assert(_graphs->insert(std::pair<ulong, Mdfg*>(_nextGraphId, newGraph)).second);
#endif

        first = newGraph->getMdfi(0);
        if(!first->setInput(next, 0)){
            throw std::runtime_error("There is an error in a MDFi link.");
        }
        /**Sends the new instruction to the interpreter.*/
#ifdef COMPUTE_COM_TIME
        unsigned long t1 = ff::getusec();
#endif
        sendToWorkers(first);
#ifdef COMPUTE_COM_TIME
        _acc += ff::getusec() - t1;
#endif
        ++_taskSent;
        ++_nextGraphId;
        ++_graphsInside;

    }
public:
    Scheduler(Mdfg *graph, InputStream *i, OutputStream *o, int parDegree,
              unsigned long int groupSize, bool orderedTasks, ff::uSWSR_Ptr_Buffer** buffers):
                _in(i),_out(o),_graph(graph),_compiled(false),_nextGraphId(0),
                _pool(NULL), _taskSent(0), _lastSent(0), _maxWorkers(parDegree),
                _numWorkers(parDegree), _groupSize(groupSize), _orderedTasks(orderedTasks),
                _buffers(buffers), _lastRcvId(0), _tasksInside(0), _graphsInside(0){
        _graphs = new std::map<ulong, Mdfg*>;
        _result = new std::map<ulong, StreamElem*>;
#ifdef COMPUTE_COM_TIME
        _acc = 0;
#endif
    }

    ~Scheduler(){
#ifdef POOL
        Mdfg* g;
        while(_pool->size()!=0){
            g = _pool->front();
            _pool->pop_front();
            delete g;
        }
        delete _pool;
#endif
        delete _graphs;
        delete _result;
    }

    void notifyRethreading(size_t oldNumWorkers, size_t newNumWorkers){
        _numWorkers = newNumWorkers;
    }


    Mdfi* schedule(){
        StreamElem* next;
        bool end = false; ///<End becomes \e true when the manager has computed all the tasks.
#ifdef POOL
        _pool = new std::deque<Mdfg*>;
        for(size_t i = 0; i < MAXPOOLSIZE; i++){
            _pool->push_back(new Mdfg(*_graph, 0));
        }
#endif

        while(!end){

            /*
            if(_graphsInside){
                std::cout << _graphsInside << " " << _tasksInside << " " << _tasksInside / _graphsInside << std::endl;
            }else{
                std::cout << _graphsInside << " " << _tasksInside << std::endl;
            }
            */

            /**Send the instructions to the interpreter.**/
            /////////////////////
            // Get from input. //
            /////////////////////
            while(keepInsertingGraphs() && _in->hasNext() && (next = _in->next()) != NULL){
                getFromInput(next);
            }

            //////////////////////
            // Get from output. //
            //////////////////////
#ifdef COMPUTE_COM_TIME
            int t1 = updateCompleted();
            _acc += t1;
#else
            updateCompleted();
#endif
            /**
             * If has received the EndOfStream and all the results
             * are calculated, then stop the manager.
             */
            if(_taskSent == 0 && !_in->hasNext()){
                end = true;
            }
            ////////////////////
            // Output stream. //
            ////////////////////
            /**
             * While exist an entry of the map with key equals to \e lastSent,
             * sends the result to the stream.
             * PRESERVES THE ORDER OF THE TASKS.
             **/
            StreamElem* se;
            std::map<ulong, StreamElem*>::iterator it;
            while((it = _result->find(_lastSent)) != _result->end()){
                se = it->second;
                _result->erase(it);
                _out->put(se);
                ++_lastSent;
            }
        }/**End of while(!end).**/
       return NULL;
#ifdef COMPUTE_COM_TIME
   std::cout << "Communication time: " << _acc << std::endl;
#endif
    }
};

WorkerMdf::WorkerMdf(ff::uSWSR_Ptr_Buffer* buffer):_buffer(buffer){;}

void WorkerMdf::compute(Mdfi* t){
    t->compute();
    _buffer->push(t);
}


Interpreter::Interpreter(Computable* c, InputStream *i, OutputStream *o, size_t parDegree,
                 unsigned long int groupSize, bool orderedTasks):
            Interpreter(compile(c), i, o, parDegree, groupSize, orderedTasks){
    //TODO Set compiledgraph
}

Interpreter::Interpreter(Mdfg *graph, InputStream *i, OutputStream *o, size_t parDegree,
                 unsigned long int groupSize, bool orderedTasks):
        _compiledGraph(NULL), _maxWorkers(parDegree){

    _p = new nornir::Parameters("parameters.xml", "archdata.xml");
    _o = new nornir::Observer;
    _p->observer = _o;

    /**Creates the SPSC queues.**/
    _buffers = new ff::uSWSR_Ptr_Buffer*[parDegree];
    for(size_t i = 0; i < parDegree; i++){
        _buffers[i] = new ff::uSWSR_Ptr_Buffer(1024);
        _buffers[i]->init();
    }

    _s = new Scheduler(graph, i, o, parDegree, groupSize, orderedTasks, _buffers);
    _farm = new nornir::Farm<Mdfi>(_p);
    _farm->addScheduler(_s);
    /**Adds the workers to the farm.**/
    for(size_t i = 0; i < parDegree; ++i){
        _farm->addWorker(new WorkerMdf(_buffers[i]));
    }
}

Interpreter::~Interpreter(){
     for(size_t i = 0; i < _maxWorkers; i++){
         delete _buffers[i];
     }
     delete[] _buffers;
     delete _p;
     delete _o;
     if(_compiledGraph)
         delete _compiledGraph;
     delete _s;
     delete _farm;
}

}
}

