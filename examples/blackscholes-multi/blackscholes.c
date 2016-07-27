// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
// 
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

// Multi-threaded pthreads header
#ifdef ENABLE_THREADS
// Add the following line so that icc 9.0 is compatible with pthread lib.
#define __thread __threadp
MAIN_ENV
#undef __thread
#endif

// Multi-threaded OpenMP header
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"

using namespace std;
using namespace tbb;
#endif //ENABLE_TBB

#ifdef ENABLE_FF
//#define BLOCKING_MODE
#include <iostream>
#include "../../src/interface.hpp"
#include "../../src/manager-multi.hpp"
#include "../../src/external/Mammut/mammut/utils.hpp"

using namespace ff;
#endif //ENABLE_FF

// Multi-threaded header for Windows
#ifdef WIN32
#pragma warning(disable : 4305)
#pragma warning(disable : 4244)
#include <windows.h>
#endif

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 100

typedef struct OptionData_ {
        fptype s;          // spot price
        fptype strike;     // strike[tid] price
        fptype r;          // risk-free interest rate[tid]
        fptype divq;       // dividend rate[tid]
        fptype v;          // volatility[tid]
        fptype t;          // time to maturity or option expiration in years 
                           //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)  
        char OptionType;   // Option type.  "P"=PUT, "C"=CALL
        fptype divs;       // dividend vals (not used in this test)
        fptype DGrefval;   // DerivaGem Reference Value
} OptionData;

OptionData **data;
fptype **prices;
int *numOptions;

int    ** otype;
fptype ** sptprice;
fptype ** strike;
fptype ** rate;
fptype ** volatility;
fptype ** otime;
int *numError;
int *nThreads;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

fptype CNDF ( fptype InputX ) 
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    } else 
        sign = 0;

    xInput = InputX;
 
    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;
    
    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal   = xLocal_1 * xNPrimeofX;
    xLocal   = 1.0 - xLocal;

    OutputX  = xLocal;
    
    if (sign) {
        OutputX = 1.0 - OutputX;
    }
    
    return OutputX;
} 

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet, int tid )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    //fptype xStockPrice;
    //fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1; 
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;    
    
    //xStockPrice = sptprice[tid];
    //xStrikePrice = strike[tid];
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log( sptprice / strike );
        
    xLogTerm = logValues;
        
    
    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;
        
    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;
    
    NofXd1 = CNDF( d1 );
    NofXd2 = CNDF( d2 );

    FutureValueX = strike * ( exp( -(rate)*(time) ) );
    if (otype == 0) {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else { 
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }
    
    return OptionPrice;
}

#ifdef ENABLE_TBB
struct mainWork {
  mainWork() {}
  mainWork(mainWork &w, tbb::split) {}

  void operator()(const tbb::blocked_range<int> &range) const {
    fptype price;
    int begin = range.begin();
    int end = range.end();

    for (int i=begin; i!=end; i++) {
      /* Calling main function to calculate option value based on 
       * Black & Scholes's equation.
       */

      price = BlkSchlsEqEuroNoDiv( sptprice[tid][i], strike[tid][i],
                                   rate[tid][i], volatility[tid][i], otime[tid][i],
                                   otype[tid][i], 0);
      prices[tid][i] = price;

#ifdef ERR_CHK 
      fptype priceDelta = data[tid][i].DGrefval - price;
      if( fabs(priceDelta) >= 1e-5 ){
        fprintf(stderr,"Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
               i, price, data[tid][i].DGrefval, priceDelta);
        numError[tid] ++;
      }
#endif
    }
  }
};

#endif // ENABLE_TBB

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_TBB
int bs_thread(void *tid_ptr) {
    int j;
    tbb::affinity_partitioner a;

    mainWork doall;
    for (j=0; j<NUM_RUNS; j++) {
      tbb::parallel_for(tbb::blocked_range<int>(0, numOptions[tid]), doall, a);
    }

    return 0;
}
#else // !ENABLE_TBB

#ifdef ENABLE_FF

// [start, end[
typedef struct{
    int start;
    int numElems;
}fftask_t;

#define CHUNKSIZE 1000

class Emitter: public nornir::Scheduler<fftask_t> {
private:
    int currentIteration;
    int currentOption;
    size_t tid;
public:
    Emitter(size_t tid):currentIteration(0), currentOption(0), tid(tid){
        ;
    }

    fftask_t* schedule(){
        while(true){
            if(currentIteration >= NUM_RUNS){
                printf("Generating end of stream.");
                fflush(stdout);
                return NULL;
            }

            fftask_t* outTask = new fftask_t;
            outTask->start = currentOption;
            outTask->numElems = 0;
            int toAssign = CHUNKSIZE;
            int availableOnIteration = numOptions[tid] - outTask->start;
            while(toAssign > 0 && currentIteration < NUM_RUNS){
                if(toAssign < availableOnIteration){
                    outTask->numElems += toAssign;
                    toAssign = 0;
                }else{
                    outTask->numElems += (availableOnIteration);
                    toAssign -= availableOnIteration;
                    availableOnIteration = numOptions[tid];
                    ++currentIteration;
                }
            }

            currentOption = (currentOption + outTask->numElems) % numOptions[tid];

            send(outTask);
        }	
    }
};

class Worker: public nornir::Worker<fftask_t> {
private:
    size_t _tid;
public:
    Worker(size_t tid):_tid(tid){;}
    void compute(fftask_t* t){
        for(uint i = 0; i < (uint) t->numElems; i++){
            fptype price;
            uint id = (t->start + i) % numOptions[_tid];
            /* Calling main function to calculate option value based on
             * Black & Scholes's equation.
             */
            price = BlkSchlsEqEuroNoDiv( sptprice[_tid][id], strike[_tid][id],
                                         rate[_tid][id], volatility[_tid][id], otime[_tid][id],
                                         otype[_tid][id], 0, _tid);
            prices[_tid][id] = price;

#ifdef ERR_CHK
            fptype priceDelta;
            priceDelta = data[_tid][id].DGrefval - price;
            if( fabs(priceDelta) >= 1e-4 ){
                printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                        id, price, data[_tid][id].DGrefval, priceDelta);
                numError[_tid] ++;
            }
#endif
        }

        delete t;
        return;
    }
};

#else // !ENABLE_FF

#ifdef WIN32
DWORD WINAPI bs_thread(LPVOID tid_ptr){
#else
int bs_thread(void *tid_ptr) {
#endif
    int i, j;
    fptype price;
#if defined(ENABLE_OPENMP) || defined(ERR_CHK)
    fptype priceDelta;
#endif
    int tid = *(int *)tid_ptr;
    int start = tid * (numOptions[tid] / nThreads[tid]);
    int end = start + (numOptions[tid] / nThreads[tid]);

    for (j=0; j<NUM_RUNS; j++) {
#ifdef ENABLE_OPENMP
#pragma omp parallel for private(i, price, priceDelta)
        for (i=0; i<numOptions[tid]; i++) {
#else  //ENABLE_OPENMP
        for (i=start; i<end; i++) {
#endif //ENABLE_OPENMP
            /* Calling main function to calculate option value based on 
             * Black & Scholes's equation.
             */
            price = BlkSchlsEqEuroNoDiv( sptprice[tid][i], strike[tid][i],
                                         rate[tid][i], volatility[tid][i], otime[tid][i],
                                         otype[tid][i], 0);
            prices[tid][i] = price;

#ifdef ERR_CHK
            priceDelta = data[tid][i].DGrefval - price;
            if( fabs(priceDelta) >= 1e-4 ){
                printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                       i, price, data[tid][i].DGrefval, priceDelta);
                numError[tid] ++;
            }
#endif
        }
    }

    return 0;
}
#endif //ENABLE_FF
#endif //ENABLE_TBB

static nornir::ManagerMulti mm;

typedef struct{
    int argc;
    char **argv;
    volatile bool started;
    size_t tid;
    nornir::ContractType ct;
    double bound;
}runArg;

void* run(void* arg){
    runArg* rarg = (runArg*) arg;
    int argc = rarg->argc;
    char** argv = rarg->argv;
    FILE *file;
    int i;
    int loopnum;
    fptype * buffer;
    int * buffer2;
    int rv;

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
        printf("PARSEC Benchmark Suite Version " __PARSEC_XSTRING(PARSEC_VERSION) "\n");
	fflush(NULL);
#else
        printf("PARSEC Benchmark Suite\n");
	fflush(NULL);
#endif //PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
   __parsec_bench_begin(__parsec_blackscholes);
#endif

   if (argc != 4)
        {
                printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
                exit(1);
        }
    nThreads[rarg->tid] = atoi(argv[1]);
    char *inputFile = argv[2];
    char *outputFile = argv[3];

    //Read input data[tid] from file
    file = fopen(inputFile, "r");
    if(file == NULL) {
      printf("ERROR: Unable to open file `%s'.\n", inputFile);
      exit(1);
    }
    rv = fscanf(file, "%i", &numOptions[rarg->tid]);
    if(rv != 1) {
      printf("ERROR: Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }
    if(nThreads[rarg->tid] > numOptions[rarg->tid]) {
      printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
      nThreads[rarg->tid] = numOptions[rarg->tid];
    }

#if !defined(ENABLE_THREADS) && !defined(ENABLE_OPENMP) && !defined(ENABLE_TBB) && !defined(ENABLE_FF)
    if(nThreads[rarg->tid] != 1) {
        printf("Error: <nthreads> must be 1 (serial version)\n");
        exit(1);
    }
#endif

    // alloc spaces for the option data[tid]
    data[rarg->tid] = (OptionData*)malloc(numOptions[rarg->tid]*sizeof(OptionData));
    prices[rarg->tid] = (fptype*)malloc(numOptions[rarg->tid]*sizeof(fptype));
    for ( loopnum = 0; loopnum < numOptions[rarg->tid]; ++ loopnum )
    {
        rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[rarg->tid][loopnum].s, &data[rarg->tid][loopnum].strike, &data[rarg->tid][loopnum].r, &data[rarg->tid][loopnum].divq, &data[rarg->tid][loopnum].v, &data[rarg->tid][loopnum].t, &data[rarg->tid][loopnum].OptionType, &data[rarg->tid][loopnum].divs, &data[rarg->tid][loopnum].DGrefval);
        if(rv != 9) {
          printf("ERROR: Unable to read from file `%s'.\n", inputFile);
          fclose(file);
          exit(1);
        }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR: Unable to close file `%s'.\n", inputFile);
      exit(1);
    }

#ifdef ENABLE_THREADS
    MAIN_INITENV(,8000000,nThreads[tid]);
#endif
    printf("Num of Options: %d\n", numOptions[rarg->tid]);
    printf("Num of Runs: %d\n", NUM_RUNS);

#define PAD 256
#define LINESIZE 64

    buffer = (fptype *) malloc(5 * numOptions[rarg->tid] * sizeof(fptype) + PAD);
    sptprice[rarg->tid] = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike[rarg->tid] = sptprice[rarg->tid] + numOptions[rarg->tid];
    rate[rarg->tid] = strike[rarg->tid] + numOptions[rarg->tid];
    volatility[rarg->tid] = rate[rarg->tid] + numOptions[rarg->tid];
    otime[rarg->tid] = volatility[rarg->tid] + numOptions[rarg->tid];

    buffer2 = (int *) malloc(numOptions[rarg->tid] * sizeof(fptype) + PAD);
    otype[rarg->tid] = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    for (i=0; i<numOptions[rarg->tid]; i++) {
        otype[rarg->tid][i]      = (data[rarg->tid][i].OptionType == 'P') ? 1 : 0;
        sptprice[rarg->tid][i]   = data[rarg->tid][i].s;
        strike[rarg->tid][i]     = data[rarg->tid][i].strike;
        rate[rarg->tid][i]       = data[rarg->tid][i].r;
        volatility[rarg->tid][i] = data[rarg->tid][i].v;
        otime[rarg->tid][i]      = data[rarg->tid][i].t;
    }

    printf("Size of data[tid]: %lu\n", numOptions[rarg->tid] * (sizeof(OptionData) + sizeof(int)));

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_begin();
#endif

#ifdef ENABLE_THREADS
#ifdef WIN32
    HANDLE *threads;
    int *nums;
    threads = (HANDLE *) malloc (nThreads[rarg->tid] * sizeof(HANDLE));
    nums = (int *) malloc (nThreads[rarg->tid] * sizeof(int));

    for(i=0; i<nThreads[rarg->tid]; i++) {
        nums[i] = i;
        threads[i] = CreateThread(0, 0, bs_thread, &nums[i], 0, 0);
    }
    WaitForMultipleObjects(nThreads[rarg->tid], threads, TRUE, INFINITE);
    free(threads);
    free(nums);
#else
    int *tids;
    tids = (int *) malloc (nThreads[tid] * sizeof(int));

    for(i=0; i<nThreads[tid]; i++) {
        tids[i]=i;
        CREATE_WITH_ARG(bs_thread, &tids[i]);
    }
    WAIT_FOR_END(nThreads[tid]);
    free(tids);
#endif //WIN32
#else //ENABLE_THREADS
#ifdef ENABLE_OPENMP
    {
        int tid=0;
        omp_set_num_threads(nThreads[tid]);
        bs_thread(&tid);
    }
#else //ENABLE_OPENMP
#ifdef ENABLE_TBB
    tbb::task_scheduler_init init(nThreads[tid]);

    int tid=0;
    bs_thread(&tid);
#else //ENABLE_TBB
#ifdef ENABLE_FF
    std::vector<ff_node*> W;
    for(int i = 0; i < nThreads[rarg->tid]; i++){
        W.push_back((ff_node*)new Worker(rarg->tid));
    }

    ff_farm<> farm(W);
    farm.add_emitter((ff_node*)new Emitter(rarg->tid));
    farm.remove_collector();

    nornir::Observer obs(std::string("stats_") + mammut::utils::intToString(rarg->tid) + std::string(".csv"), 
                         std::string("calibration_") + mammut::utils::intToString(rarg->tid) + std::string(".csv"),
                         std::string("summary_") + mammut::utils::intToString(rarg->tid) + std::string(".csv"));
    nornir::Parameters ap("parameters.xml");
    ap.contractType = rarg->ct;
    if(rarg->ct == nornir::CONTRACT_PERF_BANDWIDTH){
        ap.requiredBandwidth = rarg->bound;
    }else if(rarg->ct == nornir::CONTRACT_POWER_BUDGET){
        ap.powerBudget = rarg->bound;
    }else{
        throw std::runtime_error("Contract type not supported yet.");
    }
    ap.observer = &obs;
    ap.expectedTasksNumber = numOptions[rarg->tid] * NUM_RUNS / CHUNKSIZE;
    nornir::ManagerFarm<> amf(&farm, ap);
    mm.addManager(&amf);
    rarg->started = true;
    sleep(2);
    while(amf.running()){sleep(1);}
    std::cout << "Manager terminated." << std::endl;
    //farm.ffStats(std::cout);
#else //ENABLE_FF
    //serial version
    int ttid=0;
    bs_thread(&ttid);
#endif //ENABLE_FF
#endif //ENABLE_TBB
#endif //ENABLE_OPENMP
#endif //ENABLE_THREADS

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_end();
#endif

    //Write prices[tid] to output file
    file = fopen(outputFile, "w");
    if(file == NULL) {
      printf("ERROR: Unable to open file `%s'.\n", outputFile);
      exit(1);
    }
    rv = fprintf(file, "%i\n", numOptions[rarg->tid]);
    if(rv < 0) {
      printf("ERROR: Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }
    for(i=0; i<numOptions[rarg->tid]; i++) {
      rv = fprintf(file, "%.18f\n", prices[rarg->tid][i]);
      if(rv < 0) {
        printf("ERROR: Unable to write to file `%s'.\n", outputFile);
        fclose(file);
        exit(1);
      }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR: Unable to close file `%s'.\n", outputFile);
      exit(1);
    }

#ifdef ERR_CHK
    printf("Num Errors: %d\n", numError[tid]);
#endif
    free(data[rarg->tid]);
    free(prices[rarg->tid]);

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_bench_end();
#endif
    return NULL;
}


#define NUM_INSTANCES 2
int main(int argc, char **argv){
    data = (OptionData**) malloc(sizeof(OptionData*)*NUM_INSTANCES);
    prices = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    numOptions = (int*) malloc(sizeof(int)*NUM_INSTANCES);

    otype = (int**) malloc(sizeof(int*)*NUM_INSTANCES);
    sptprice = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    strike = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    rate = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    volatility = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    otime = (fptype**) malloc(sizeof(fptype*)*NUM_INSTANCES);
    numError = (int*) malloc(sizeof(int)*NUM_INSTANCES);
    for(size_t i = 0; i < NUM_INSTANCES; i++){
        numError[i] = 0;
    }
    nThreads = (int*) malloc(sizeof(int)*NUM_INSTANCES);

    size_t nexttid = 0;
    pthread_t tid1, tid2;
    mm.start();
    runArg* rarg = new runArg();
    rarg->argc = argc;
    rarg->argv = argv;
    rarg->started = false;
    rarg->tid = nexttid;
    rarg->ct = nornir::CONTRACT_PERF_BANDWIDTH;
    rarg->bound = 14420;
    nexttid++;
    std::cout << "Creating first blackscholes instance." << std::endl;
    pthread_create(&tid1, NULL, &run, (void*) rarg);
    std::cout << "Instance created." << std::endl;
    while(!rarg->started){;}
    std::cout << "Instance started." << std::endl;
    sleep(8);

    rarg = new runArg();
    rarg->argc = argc;
    rarg->argv = argv;
    rarg->started = false;
    rarg->tid = nexttid;
    rarg->ct = nornir::CONTRACT_POWER_BUDGET;
    rarg->bound = 50;
    //    rarg->bound = 14420*2;
    nexttid++;
    std::cout << "Creating second blackscholes instance." << std::endl;
    pthread_create(&tid2, NULL, &run, (void*) rarg);
    std::cout << "Instance created." <<std::endl;
    while(!rarg->started){;}
    std::cout << "Instance started." <<std::endl;

    pthread_join(tid1, NULL);
    pthread_join(tid2, NULL);

    return 0;
}

