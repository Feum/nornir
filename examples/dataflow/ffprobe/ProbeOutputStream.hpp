/*
 * ProbeOutputStream.hpp
 *
 *  \date 14/mag/2010
 *  \date 09/ago/2010
 *  \author Daniele De Sensi (d.desensi@gmail.com)
 */

#ifndef PROBEOUTPUTSTREAM_HPP_
#define PROBEOUTPUTSTREAM_HPP_

#include "../../../src/dataflow/stream.hpp"

namespace faskelProbe{
class ProbeOutputStream: public nornir::dataflow::OutputStream{
private:
    FILE* out; ///<File where to print the flows in textual format.
    uint qTimeout, ///<It specifies how long expired flows (queued before delivery) are emitted
        flowSequence, ///<Sequence number for the flows to export
        minFlowSize;///<Minimum tcp flows size
    std::queue<hashElement*>* q; ///< Queue of expired flows
    time_t lastEmission; ///< Time of the last export
    Exporter ex;
    bool _deallocTasks;
    /**
     * Exports the flow to the remote collector (also prints it into the file).
     */
    inline void exportFlows(){
        int size=q->size();
        ex.sendToCollector(q,flowSequence,out);
        flowSequence+=size;
    }
public:
    /**
     * Constructor of the output stream.
     * \param out The FILE* where to print exported flows.
     * \param queueTimeout It specifies how long expired flows (queued before delivery) are emitted.
     * \param collector The host of the collector.
     * \param port The port where to send the flows.
     * \param minFlowSize If a TCP flow doesn't have more than minFlowSize bytes isn't exported (0 is unlimited).
     * \param systemStartTime The system start time.
     * \param deallocTasks If true, tasks will be automatically deallocated.
     */
    inline ProbeOutputStream(FILE* out, uint queueTimeout, char const* collector,
                             uint port, uint minFlowSize, uint32_t systemStartTime,
                             bool deallocTasks = true):
            out(out),qTimeout(queueTimeout),flowSequence(0),
            minFlowSize(minFlowSize),q(new std::queue<hashElement*>),lastEmission(time(NULL)),
            ex(collector,port,ffalloc,systemStartTime), _deallocTasks(deallocTasks){
        if(out!=NULL)
            fprintf(out,"IPV4_SRC_ADDR|IPV4_DST_ADDR|OUT_PKTS|OUT_BYTES|FIRST_SWITCHED|LAST_SWITCHED|L4_SRC_PORT|L4_DST_PORT|TCP_FLAGS|"
                    "PROTOCOL|SRC_TOS|\n");
        struct sigaction s;
        bzero( &s, sizeof(s) );
        s.sa_handler=handler;
        sigaction(SIGINT,&s,NULL);
        assert(ffalloc!=NULL);
        if (ffalloc->register4free()<0) {
            fprintf(stderr,"lastStage, register4free fails\n");
            exit(-1);
        }
    }

    /**
     * Destructor of the output stream.
     */
    inline ~ProbeOutputStream(){
        delete q;
    }

    /**
     * Puts the task into the output stream.
     * \param a The task to put.
     */
    inline void put(void* a){
        if(!_deallocTasks){
            return;
        }
        ProbeTask* t=(ProbeTask*) a;
        hashElement* f;
        if(a!=NULL){
            myList<hashElement*>* l = t->getFlowsToExport();
            assert(t->elementsToAddSize() == 0);
            time_t now = time(NULL);
            while(l->size()!=0){
                assert(l->pop(&f)==0);
                if(!(f->prot==TCP_PROT_NUM && f->dOctets<minFlowSize))
                    q->push(f);
                if(q->size()==30){
                    exportFlows();
                    lastEmission = now;
                }
            }
            /**Exports flows every qTimeout seconds.**/
            if((t->isEof() || (now - lastEmission>=qTimeout)) && !q->empty()){
                exportFlows();
                lastEmission = now;
            }
            if(_deallocTasks){
                ffalloc->free(t);
            }
        }else{
            if(!q->empty())
                exportFlows();
        }
    }
};
}
#endif /* PROBEOUTPUTSTREAM_HPP_ */
