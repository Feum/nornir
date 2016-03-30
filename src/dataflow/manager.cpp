/*
 * manager.cpp
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

#include "manager.hpp"
#include "../external/Mammut/mammut/mammut.hpp"

using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

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
        Pipeline* p=(Pipeline*) c;
        /**Compiles the two stages.**/
        Mdfg* g1=compile(p->getFirstStage());
        Mdfg* g2=compile(p->getSecondStage());
        bool lastSet=false;
        int n=g2->getNumMdfi();
        Mdfi *added,*toInsert,
            *oldSecondStageFirst,///<The first instruction of the second stage of the pipeline.
            *oldLast,///<The last instruction of the first stage of the pipeline.
            *last;///<The last instruction of the new graph (Stage1 & Stage2).

        /**
         * \e newIds is an array of Mdfi identifiers. \e newIds[i] is the new
         * identifier of the instruction of the second stage, that previously
         * has \e i as identifier.
         **/
        int *newIds=new int[n];
        oldLast=g1->getLast();
        oldSecondStageFirst=g2->getFirst();
        /**
         * Discriminates the case where the second stage has only one
         * instruction.
         **/
        if(n==1){
            oldSecondStageFirst=g1->createLastMdfi(oldSecondStageFirst);
            lastSet=true;
        }else
            oldSecondStageFirst=g1->createMdfi(oldSecondStageFirst);
        newIds[0]=oldSecondStageFirst->getId();
        /**
         * Adds the instructions of the second stage (except first and last
         * instruction).
         **/
        for(int i=1; i<n-1;i++){
            toInsert=g2->getMdfi(i);
            added=g1->createMdfi(toInsert);
            newIds[i]=added->getId();
        }

        /**If the last instruction isn't set.**/
        if(!lastSet){
            last=g2->getLast();
            last=g1->createLastMdfi(last);
            newIds[n-1]=last->getId();
        }
        /**Updates destinations.**/
        for(int i=0; i<n; i++)
            g1->getMdfi(newIds[i])->updateDestinations(newIds);

        /**Links the two stages.**/
        g1->link(oldLast,0,oldSecondStageFirst,0);

        delete[] newIds;
        /**
         * Deletes the old graph of the second stage because all its
         * instructions were copied in the graph of the first stage.
         **/
        delete g2;
        return g1;
    /**Map and reduce compiling.**/
    }else if(tc==typeid(EmitterWorkerCollector)){
        EmitterWorkerCollector* ewc=(EmitterWorkerCollector*) c;
        int workersNum=ewc->getNWorkers();
        /**Compiles worker.**/
        Mdfg* worker=compile(ewc->getWorker());
        Mdfg* emitter=new Mdfg(ewc->getEmitter(),1,workersNum);
        Mdfg* collector=new Mdfg(ewc->getCollector(),workersNum,1);

        int workerSize=worker->getNumMdfi();
        /**
         * \e newWorkerIds is an array of Mdfi identifiers. \e newWorkerIds[i]
         * is the new identifier of the instruction of the worker, that
         * previously has \e i as identifier.
         **/
        int *newWorkerIds=new int [workerSize];
        /**
         * \e firstWorkerInstr and \e lastWorkerInstrs contain the first and
         * the last instructions of the workers.
         **/
        int *firstWorkerInstr=new int[workersNum],*lastWorkerInstr=new int[workersNum];
        /**Adds the workers.**/
        for(int i=0; i<workersNum; i++){
            for(int j=0; j<workerSize; j++)
                newWorkerIds[j]=emitter->createMdfi(worker->getMdfi(j))->getId();
            firstWorkerInstr[i]=newWorkerIds[0];
            lastWorkerInstr[i]=newWorkerIds[workerSize-1];
            /**Updates the destinations.**/
            for(int j=0; j<workerSize; j++)
                emitter->getMdfi(newWorkerIds[j])->updateDestinations(newWorkerIds);
        }
        delete[] newWorkerIds;
        /**Links the emitter to the workers.**/
        for(int i=0; i<workersNum; i++)
            emitter->link(emitter->getFirst(),i,emitter->getMdfi(firstWorkerInstr[i]),0);
        delete[] firstWorkerInstr;
        /**Adds the collector.**/
        Mdfi *added=collector->getFirst(),*firstCollInstr;
        firstCollInstr=emitter->createLastMdfi(added);

        /**Links the workers to collector.**/
        for(int i=0; i<workersNum; i++)
            emitter->link(emitter->getMdfi(lastWorkerInstr[i]),0,firstCollInstr,i);
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

Manager::Manager(Computable* c, InputStream *i, OutputStream *o, int parDegree,
                 unsigned long int groupSize, bool orderedTasks):
            in(i),out(o),compiled(true),nextGraphId(0),groupSize(groupSize),
            executed(0),taskSent(0),lastSent(0){
    graph=compile(c);
    intr=new Interpreter(parDegree, orderedTasks);

#ifdef POOL
    pool=new std::deque<Mdfg*>;
#endif
    fireable=new std::deque<Mdfi*>;
    graphs=new hashMap<Mdfg*>;
    result=new hashMap<StreamElem*>;
    tempTask=new StreamElem*;
#ifdef COMPUTE_COM_TIME
    acc=0;
#endif
}

Manager::Manager(Mdfg *graph, InputStream *i, OutputStream *o, int parDegree,
                 unsigned long int groupSize, bool orderedTasks):
            in(i),out(o),graph(graph),compiled(false),nextGraphId(0),
            groupSize(groupSize),executed(0),taskSent(0),lastSent(0){
    intr=new Interpreter(parDegree, orderedTasks);
#ifdef POOL
    pool=new std::deque<Mdfg*>;
#endif
    fireable=new std::deque<Mdfi*>;
    graphs=new hashMap<Mdfg*>;
    result=new hashMap<StreamElem*>;
    tempTask=new StreamElem*;
#ifdef COMPUTE_COM_TIME
    acc=0;
#endif
}

Manager::~Manager(){
    delete intr;
#ifdef POOL
    Mdfg* g;
    while(pool->size()!=0){
        g=pool->front();
        pool->pop_front();
        delete g;
    }
    delete pool;
#endif
    delete graphs;
    delete result;
    delete fireable;
    if(compiled)
        delete graph;
    delete tempTask;
}

void Manager::stats(std::ostream& out){
    intr->stats(out);
}

void Manager::exec(){
    StreamElem* temp;
    bool end=false; ///<End becomes \e true when the manager has computed all the tasks.
    ff::squeue<OutputToken> res; ///<Queue of results computed by the interpreter.
    Mdfg **g=new Mdfg*;
    Mdfg **tdl=new Mdfg*;
    Mdfi *ins;
    OutputToken ot;
    TokenId dest;
    unsigned long int graphId;

    Mammut m;
    TasksManager* pm = m.getInstanceTask();

    ProcessHandler* process = pm->getProcessHandler(getpid());
    ThreadHandler* thread = process->getThreadHandler(gettid());
    thread->move((VirtualCoreId) 23);
    process->releaseThreadHandler(thread);
    pm->releaseProcessHandler(process);

    while(!end){
        /**Send the instructions to the interpreter.**/
        while(executed<groupSize && fireable->size()){
            intr->exec(fireable->front());
            fireable->pop_front();
            ++executed;
        }
        /**Check if there are new tasks on the input stream.**/
        getFromInput();
        /**Awaits the results.**/
#ifdef COMPUTE_COM_TIME
        int t1=intr->wait(res);
        acc+=t1;
#else
        intr->wait(res);
#endif
        executed=0;
        /**Updates the destinations.**/
        while(res.size()){
            ot=res.back();
            res.pop_back();
            dest=ot.getDest();
            graphId=dest.getGraphId();
            /**
             * If was the last instruction of a graph's copy, puts it
             * into the output vector and delete
             * the copy of the graph.
             */
            if(dest.isOutStream()){
                result->put(graphId,ot.getResult());
                graphs->get(graphId,tdl);
                graphs->erase(graphId);
#ifdef POOL
                if(pool->size()<MAXPOOLSIZE)
                    pool->push_back(*tdl);
                else
                    delete *tdl;
#else
                delete *tdl;
#endif
                --taskSent;
                /**
                 * If has received the EndOfStream and all the results
                 * are calculated, then stop the manager.
                 */
                if(taskSent == 0 && !in->hasNext())
                    end=true;
            }else{
                /**Takes the pointer to the copy of the graph.**/
                graphs->get(graphId,g);
                /**Takes the pointer to the instruction.**/
                ins=(*g)->getMdfi(dest.getMdfId());
                /**Updates the instruction adding the input token.**/
                temp=ot.getResult();
                ins->setInput(temp,dest.getTokId());
                /**
                 * If the instruction is fireable, adds it to the pool of
                 * fireable instructions.
                 **/
                if(ins->isFireable()){
                    if(executed<groupSize){
                        intr->exec(ins);
                        ++executed;
                    }else
                        fireable->push_back(ins);
                }
            }
        }
        /**Sends new results to the output stream.**/
        flushOnStream();
    }/**End of while(!end).**/
       intr->stop();
       delete g;
       delete tdl;
#ifdef COMPUTE_COM_TIME
       std::cout << "Communication time: " << acc << std::endl;
#endif
}

void Manager::getFromInput(){
    Mdfg *newGraph;
    Mdfi *first;
    StreamElem* next;
    /**
     * Tries to read from the stream while the number of task submitted to
     * the interpreter is less than \e groupSize
     * and while there are other task on the input stream.
     */

    while(in->hasNext() &&  executed < groupSize && (next=in->next())!=NULL ){
#ifdef POOL
        if(pool->size()!=0){
            newGraph=pool->front();
            pool->pop_front();
            newGraph->reset(nextGraphId);
        }else
            newGraph=new Mdfg(*graph,nextGraphId);
#else
        newGraph=new Mdfg(*graph,nextGraphId);
#endif
        graphs->put(nextGraphId,newGraph);
        first=newGraph->getFirst();
        if(!first->setInput(next,0)){
            std::cerr << "There is an error in a MDFi link." << std::endl;
            exit(-1);
        }
        /**Sends the new instruction to the interpreter.*/
#ifdef COMPUTE_COM_TIME
        unsigned long t1=ff::getusec();
#endif
        intr->exec(first);
#ifdef COMPUTE_COM_TIME
        acc+=ff::getusec()-t1;
#endif
        ++taskSent;
        ++executed;
        ++nextGraphId;
    }
}

void Manager::flushOnStream(){
    /**
     * While exist an entry of the map with key equals to \e lastSent,
     * sends the result to the stream.
     **/
    while(result->getAndErase(lastSent,tempTask)){
        out->put(*tempTask);
        ++lastSent;
    }
}

}
}
