/*
 * ProbeInputStream.hpp
 *
 *  \date 14/mag/2010
 *  \date 09/ago/2010
 *  \author Daniele De Sensi (d.desensi@gmail.com)
 */

#ifndef PROBEINPUTSTREAM_HPP_
#define PROBEINPUTSTREAM_HPP_

#include "../../../src/dataflow/stream.hpp"


namespace faskelProbe{

class ProbeInputStream: public nornir::dataflow::InputStream{
private:
    int cnt, ///< Maximum number of packet to read from the device (or from the .pcap file)
        pktRcvd; ///< Packet received after a call of 'pcap_dispatch'
    bool offline, ///< True if the device is a .pcap file
        end; ///< Flag for the termination of the pipeline
public:
    /**
     * Constructor of the input stream.
     * \param nw Number of workers.
     * \param device Name of the device (or of the .pcap file)
     * \param noPromisc False if the interface must be set in promiscous mode, false otherwise.
     * \param filter_exp A bpf filter.
     * \param cnt Maximum number of packet to read from the device (or from the .pcap file)
     * \param h Size of the hash table (Sizeof(HashOfWorker1)+Sizeof(HashOfWorker2)+...+Sizeof(HashOfWorkerN))
     * \param readTimeout The read timeout when reading from the pcap socket
     * \param alloc A pointer to the fastflow's allocator
     */
    inline ProbeInputStream(int nw, const char* device, bool noPromisc, char* filter_exp, int cnt, int h, int readTimeout, ff::ff_allocator *alloc):
    cnt(cnt),pktRcvd(0),end(false){
        ffalloc=alloc;
        nWorkers=nw!=0?nw:1;
        quit=false;
        hsize=h;
        char errbuf[PCAP_ERRBUF_SIZE];
        struct bpf_program fp;      /* The compiled filter expression */
        /**Accepts only ipv4 traffic.**/
        bpf_u_int32 mask;       /* The netmask of our sniffing device */
        bpf_u_int32 net;        /* The IP of our sniffing device */
        if (pcap_lookupnet(device, &net, &mask, errbuf) == -1) {
            handle=pcap_open_offline(device,errbuf);
            offline=true;
            net = 0;
            mask = 0;
        }else{
            int prom=noPromisc?0:1;
            handle = pcap_open_live(device, 200, prom, readTimeout, errbuf);
            offline=false;
        }
        if (handle == NULL) {
            fprintf(stderr, "Couldn't open device %s: %s\n", device, errbuf);
            exit(-1);
        }
        /**Accepts only ipv4 traffic (netflow v5 doesn't support other network protocols).**/
        if (pcap_compile(handle, &fp, "ip", 0, net) == -1) {
            fprintf(stderr, "Couldn't parse filter %s: %s\n", filter_exp, pcap_geterr(handle));
            exit(-1);
        }

        if (pcap_setfilter(handle, &fp) == -1) {
            fprintf(stderr, "Couldn't install filter %s: %s\n", filter_exp, pcap_geterr(handle));
            exit(-1);
        }

        pcap_freecode(&fp);

        /**Sets the filter passed by user.**/
        if(filter_exp!=NULL){
            if (pcap_compile(handle, &fp, filter_exp, 0, net) == -1) {
                fprintf(stderr, "Couldn't parse filter %s: %s\n", filter_exp, pcap_geterr(handle));
                exit(-1);
            }
            if (pcap_setfilter(handle, &fp) == -1) {
                fprintf(stderr, "Couldn't install filter %s: %s\n", filter_exp, pcap_geterr(handle));
                exit(-1);
            }
            pcap_freecode(&fp);
        }
        int datalinkType=pcap_datalink(handle);
        //TODO Add other switch-case to add the support to other datalink's protocols.
        switch(datalinkType){
            case 1:
                datalinkOffset=14;
                break;
            default:
                fprintf(stderr, "Datalink offset for datalink type: %d unknown.",datalinkType);
                exit(-1);
        }
        /**
         * Signal handling.
         */
        struct sigaction s;
        bzero( &s, sizeof(s) );
        s.sa_handler=handler;
        sigaction(SIGINT,&s,NULL);
        /**
         * Registers this thread as allocator
         */
        if (ffalloc->registerAllocator()<0){
            fprintf(stderr,"ffalloc->registerAllocator() fail.");
            exit(-1);
        }
    }

    /**
     * Destructor of the input stream.
     */
    inline ~ProbeInputStream(){
        ffalloc->deregisterAllocator();
        pcap_close(handle);
    }

    /**
     * Returns the next element of the stream.
     * \return The next element of the stream.
     */
    nornir::dataflow::Task* next(){
        ProbeTask* t=(ProbeTask*) ffalloc->malloc(sizeof(ProbeTask));
        t->init(nWorkers,ffalloc);
        pktRcvd=pcap_dispatch(handle,cnt,dispatchCallback,(u_char*)t);
        if((pktRcvd==0 && offline) || quit){
            end=true;
            t->setEof();
        }
        else if(pktRcvd==0 && !offline){
            t->setReadTimeoutExpired();
        }
        return t;
    }

    /**
     * Checks if the EndOfStream is arrived.
     * \return \e True if the EndOfStream is arrived, \e false otherwise.
     */
    inline bool hasNext(){
        return !end;
    }
};
}
#endif /* PROBEINPUTSTREAM_HPP_ */
