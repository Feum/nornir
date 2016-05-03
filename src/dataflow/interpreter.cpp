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
#include "../external/Mammut/mammut/mammut.hpp"
#include <map>

using namespace mammut;
using namespace mammut::topology;

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
    if(tc == typeid(Farm)){
        Farm* f = (Farm*) c;
        Mdfg* g = compile(f->getWorker());
        g->init();
        return g;
    /**
     * Pipeline compiling. Adds all istructions of the second stage in the
     * first stage of the pipeline and links it together.
     **/
    }else if(tc == typeid(Pipeline)){
        Pipeline* p = (Pipeline*) c;
        /**Compiles the two stages.**/
        Mdfg* g1 = compile(p->getFirstStage());
        Mdfg* g2 = compile(p->getSecondStage());
        bool lastSet = false;
        uint n = g2->getNumMdfi();
        uint added, oldSecondStageFirstId;
        Computable *oldSecondStageFirst, ///<The first instruction of the second stage of the pipeline.
                   *oldLast; ///<The last instruction of the first stage of the pipeline.

        /**
         * \e newIds is an array of Mdfi identifiers. \e newIds[i] is the new
         * identifier of the instruction of the second stage, that previously
         * has \e i as identifier.
         **/
        int *newIds = new int[n];
        oldLast = g1->getLast();
        oldSecondStageFirst = g2->getFirst();
        oldSecondStageFirstId = g1->createMdfi(g2, g2->getFirstId());
        /**
         *  Ok to call getMdfi(...) even if the graph will me modified. Indeed
         *  we will not store the Mdfi itself but only the computable associated
         *  to it.
         */
        oldSecondStageFirst = g1->getMdfi(oldSecondStageFirstId)->getComputable();
        /**
         * Discriminates the case where the second stage has only one
         * instruction.
         **/
        if(n == 1){
            lastSet = true;
        }
        newIds[0] = oldSecondStageFirstId;
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
            newIds[n - 1] = g1->createMdfi(g2, g2->getLastId());
        }
        /**Updates destinations.**/
        for(uint i = 0; i < n; i++){
            g1->updateDestinations(newIds[i], newIds);
        }

        /**Links the two stages.**/
        g1->link(oldLast, oldSecondStageFirst);

        delete[] newIds;
        /**
         * Deletes the old graph of the second stage because all its
         * instructions were copied in the graph of the first stage.
         **/
        delete g2;

        g1->init();
        return g1;
    /**Map and reduce compiling.**/
    }else if(tc == typeid(EmitterWorkerCollector)){
        EmitterWorkerCollector* ewc = (EmitterWorkerCollector*) c;
        int workersNum = ewc->getNWorkers();
        /**Compiles worker.**/
        Mdfg* worker = compile(ewc->getWorker());
        Mdfg* emitter = new Mdfg(ewc->getEmitter());
        Mdfg* collector = new Mdfg(ewc->getCollector());

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
            emitter->link(emitter->getFirst(),
                          emitter->getMdfi(firstWorkerInstr[i])->getComputable());
        }
        delete[] firstWorkerInstr;
        /**Adds the collector.**/
        uint firstCollInstr = emitter->createMdfi(collector, 0);

        /**Links the workers to collector.**/
        for(int i = 0; i < workersNum; i++){
            emitter->link(emitter->getMdfi(lastWorkerInstr[i])->getComputable(),
                          emitter->getMdfi(firstCollInstr)->getComputable());
        }
        delete[] lastWorkerInstr;
        /**
         * Deletes the worker and the collector because all their instructions
         *  were copied in the graph of the emitter.
         **/
        delete worker;
        delete collector;

        emitter->init();
        return emitter;
    /**If c is an unknown skeleton, then it is a sequential skeleton.**/
    }else{
        Mdfg* g = new Mdfg(c);
        g->init();
        return g;
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
    std::map<ulong, void*> *_result;

    unsigned long int
    /**Number of results not yet calculated.**/
        _taskSent,
    /**Index of the last task sent to the output stream.**/
        _lastSent;
    size_t _maxWorkers;
    size_t _numWorkers;
    size_t _maxGraphs;
    bool _orderedProc;
    bool _orderedOut;
    QUEUE& _q;
    size_t _lastRcvId;
    ulong _tasksInside;
    ulong _graphsInside;
    size_t* _scheduling;
    size_t* _numMdfi;
    size_t* _mdfiSent;
    std::vector<WorkerMdf*>& _workers;

    inline void sendToWorkers(Mdfi* instr){
        size_t destination;
        ++_tasksInside;
        if(_orderedProc){
            destination = _scheduling[instr->getId()];
            ++_mdfiSent[destination];
            sendTo(instr, destination);
        }else{
            send(instr);
        }
    }

    void updateScheduling(size_t insId){
#if 0
        size_t minId = 0;
        size_t minLoad = std::numeric_limits<size_t>::max();
        size_t load;
        size_t baseId = _scheduling[insId];
        size_t workerId = 0;
        for(size_t i = 0; i < _numWorkers; i++){
            workerId = (baseId + i) % _numWorkers;
            load = _mdfiSent[workerId] - _workers.at(workerId)->getProcessedTasks(); // Approximates the number of tasks in the input Q of worker i.
            //std::cout << "Load of " << workerId << ": " << load << std::endl;
            if(load < minLoad){
                minLoad = load;
                minId = workerId;
            }
        }
        _scheduling[insId] = minId;
#endif
    }

    void getFromInput(void* next){
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

        first = newGraph->getMdfi(newGraph->getFirstId());
        first->setInput(next, NULL);
        ++_numMdfi[first->getId()];
        /**Sends the new instruction to the interpreter.*/
        sendToWorkers(first);
        ++_taskSent;
        ++_nextGraphId;
        ++_graphsInside;

    }
public:
    Scheduler(Mdfg *graph, InputStream *i, OutputStream *o, size_t parDegree,
              QUEUE& q, std::vector<WorkerMdf*>& workers, Parameters* p):
                _in(i), _out(o), _graph(graph), _compiled(false), _nextGraphId(0),
                _pool(NULL), _taskSent(0), _lastSent(0), _maxWorkers(parDegree),
                _numWorkers(parDegree), _maxGraphs(p->dataflow.maxGraphs),
                _orderedProc(p->dataflow.orderedProcessing),
                _orderedOut(p->dataflow.orderedOutput), _q(q), _lastRcvId(0),
                _tasksInside(0), _graphsInside(0), _workers(workers){
        _graphs = new std::map<ulong, Mdfg*>;
        _result = new std::map<ulong, void*>;
        _scheduling = new size_t[_graph->getNumMdfi()];
        _numMdfi = new size_t[_graph->getNumMdfi()];
        for(size_t i = 0; i < _graph->getNumMdfi(); i++){
            _scheduling[i] = i % parDegree;
            _numMdfi[i] = 0;
        }
        _mdfiSent = new size_t[parDegree];
        for(size_t i = 0; i < parDegree; i++){
            _mdfiSent[i] = 0;
        }
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
        delete[] _scheduling;
        delete[] _numMdfi;
        delete[] _mdfiSent;
    }

    void notifyRethreading(size_t oldNumWorkers, size_t newNumWorkers){
        _numWorkers = newNumWorkers;
        for(size_t i = 0; i < _graph->getNumMdfi(); i++){
            _scheduling[i] = i % _numWorkers;
            _numMdfi[i] = 0;
        }
    }


    Mdfi* schedule(){
        void* next;
        bool end = false; ///<End becomes \e true when the manager has computed all the tasks.
        OutputToken ot;
        ulong graphId;
        TokenId dest;
        Mdfi* poppedIns;
        Mdfi* ins;
        void* task;
        uint poppedId;
        _q.registerq(0);

#ifdef POOL
        _pool = new std::deque<Mdfg*>;
        for(size_t i = 0; i < MAXPOOLSIZE; i++){
            _pool->push_back(new Mdfg(*_graph, 0));
        }
#endif

        /** Bootstrap. **/
        size_t inserted = 0;
        size_t totalSent = 0;
        while(inserted < _maxGraphs && _in->hasNext()){
            next = _in->next();
            if(next){
                getFromInput(next);
                ++inserted;
                ++totalSent;
            }
        }
        uint popped = 0;
        while(!end){
            /////////////////////
            // Get from input. //
            ///////////////////// 
            if(_graphsInside < _maxGraphs &&
               _in->hasNext() && (next = _in->next()) != NULL){
                getFromInput(next);
                //popped = 0;
                ++inserted;
                ++totalSent;
            }

            //////////////////////
            // Get from output. //
            //////////////////////
            if(_q.pop((void**) &task, 0)){
                popped++;
                --_tasksInside;
                poppedIns = ((Mdfi*) task);
                for(uint j = 0; j < poppedIns->getNumOutTokens(); j++){
                    ot = *(poppedIns->getOutToken(j));
                    poppedId = poppedIns->getId();
                    --_numMdfi[poppedId];
                    if(_numMdfi[poppedId] == 0){
                        // Can change the scheduling.
                        updateScheduling(poppedId);
                    }
                    dest = ot.getDest();
                    graphId = dest.getGraphId();
                    Mdfg* currentGraph = NULL;
                    /**
                     * If was the last instruction of a graph's copy, puts it
                     * into the output vector and delete
                     * the copy of the graph.
                     */
                    if(dest.isOutStream()){
                        if(_out){
                            if(!_orderedOut || graphId == _lastSent){
                                _out->put(ot.getResult());
                                ++_lastSent;
                            }else{
                                assert(_result->emplace(std::piecewise_construct,
                                       std::forward_as_tuple(graphId),
                                       std::forward_as_tuple(ot.getResult())).second);
                            }
                        }
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
                        ins->setInput(ot.getResult(), poppedIns->getComputable());
                        /**
                         * If the instruction is fireable, adds it to the pool of
                         * fireable instructions.
                         **/
                        if(ins->isFireable()){
                            sendToWorkers(ins);
                            ++_numMdfi[ins->getId()];
                            ++totalSent;
                        }
                    }
                }
            }

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
            if(_out && _orderedOut){
                void* se;
                std::map<ulong, void*>::iterator it;
                while((it = _result->find(_lastSent)) != _result->end()){
                    se = it->second;
                    _result->erase(it);
                    _out->put(se);
                    ++_lastSent;
                }
            }
        }/**End of while(!end).**/
        _q.deregisterq(0);
        return NULL;
    }
};

WorkerMdf::WorkerMdf(QUEUE& q):_q(q), _init(false), _qId(-1), _processedTasks(0){;}

void WorkerMdf::compute(Mdfi* t){
    if(!_init){
        _qId = getId() + 1;
        _init = true;
        _q.registerq(_qId);
    }
    t->compute();
    _q.push((void*) t, _qId);
    ++_processedTasks;
}

size_t WorkerMdf::getProcessedTasks() const{
    return _processedTasks;
}


Interpreter::Interpreter(Parameters* p, Computable* c, InputStream *i, OutputStream *o):
            Interpreter(p, compile(c), i, o){
    //TODO Set compiledgraph
}

Interpreter::Interpreter(Parameters* p, Mdfg *graph, InputStream *i, OutputStream *o):
        _p(p), _compiledGraph(NULL){
    graph->init();
    Mammut m;
    size_t numPhysicalCores = m.getInstanceTopology()->getPhysicalCores().size();
    if(numPhysicalCores < 3){
        throw std::runtime_error("Not enough cores available (you need at least "
                                 "3 physical cores).");
    }
    /** -2: One for the dataflow scheduler and one for nornir manager. **/
    _maxWorkers = numPhysicalCores - 2;

    /**Creates the SPSC queues.**/
    _q.init(_maxWorkers + 1); /* +1 for the scheduler. */

    _s = new Scheduler(graph, i, o, _maxWorkers, _q, _workers, p);
    _farm = new nornir::Farm<Mdfi>(_p);
    _farm->addScheduler(_s);
    /**Adds the workers to the farm.**/
    for(size_t i = 0; i < _maxWorkers; ++i){
        _workers.push_back(new WorkerMdf(_q));
        _farm->addWorker(_workers.back());
    }
}

Interpreter::~Interpreter(){
     if(_compiledGraph){
         delete _compiledGraph;
     }
     delete _s;
     delete _farm;
     for(size_t i = 0; i < _workers.size(); i++){
         delete _workers.at(i);
     }
}

}
}

