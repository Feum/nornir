#include <cstdlib>
#include <cstdio>
#include <vector>
#include <farm.hpp>
#include <mammut/mammut.hpp>

using namespace ff;


class Obs: public adpff::adp_ff_farm_observer{
private:
  std::ofstream _statsFile;
  std::ofstream _energyFile;
  mammut::energy::JoulesCpu _totalUsedJoules, _totalUnusedJoules;
  unsigned long samples;
public:
  Obs():samples(0){
    _statsFile.open("stats.txt");
    if(!_statsFile.is_open()){
      throw std::runtime_error("Obs: Impossible to open stats file.");
    }
    _statsFile << "# [[EmitterVc][WorkersVc][CollectorVc]] NumWorkers,Frequency CurrentBandwidth CurrentUtilization" << std::endl;

    _energyFile.open("energy.txt");
    if(!_energyFile.is_open()){
      throw std::runtime_error("Obs: Impossible to open energy file.");
    }
    _energyFile << "# UsedVCCpuEnergy UsedVCCoresEnergy UsedVCGraphicEnergy UsedVCDRAMEnergy NotusedVCCpuEnergy NotusedVCCoresEnergy NotusedVCGraphicEnergy NotusedVCDRAMEnergy " << std::endl;
  }

  ~Obs(){
    _energyFile << _totalUsedJoules.cpu/(double)samples << " " << _totalUsedJoules.cores/(double)samples << " " << _totalUsedJoules.graphic/(double)samples << " " << _totalUsedJoules.dram/(double)samples << std::endl;
    _statsFile.close();
    _energyFile.close();
  }

  void observe(){
    /****************** Stats ******************/
    _statsFile << time(NULL);
    _statsFile << " [";
    if(_emitterVirtualCore){
      _statsFile << "[" << _emitterVirtualCore->getVirtualCoreId() << "]";
    }

    _statsFile << "[";
    for(size_t i = 0; i < _workersVirtualCore.size(); i++){
      _statsFile << _workersVirtualCore.at(i)->getVirtualCoreId() << ",";
    }
    _statsFile << "]";

    if(_collectorVirtualCore){
      _statsFile << "[" << _collectorVirtualCore->getVirtualCoreId() << "]";
    }
    _statsFile << "] ";

    _statsFile << _numberOfWorkers << "," << _currentFrequency << " ";
    _statsFile << _currentBandwidth << " ";
    _statsFile << _currentUtilization << " ";
    _statsFile << _currentCoresUtilization << " ";

    _statsFile << std::endl;
    /****************** Energy ******************/

    _energyFile << _usedJoules.cpu << " " << _usedJoules.cores << " " << _usedJoules.graphic << " " << _usedJoules.dram << " ";
    _energyFile << _unusedJoules.cpu << " " << _unusedJoules.cores << " " << _unusedJoules.graphic << " " << _unusedJoules.dram << " ";
    _energyFile << std::endl;
    _totalUsedJoules += _usedJoules;
    _totalUnusedJoules += _unusedJoules;
    ++samples;
  }
};



/* globals */
double* A = NULL;
double* B = NULL;
double* C = NULL;

void Print(int N, int size) {
    int tsize=size*size;
    for(int i=0;i<N;++i) {
	for(int j=0;j<size;++j) {
	    for(int k=0;k<size;++k) {
		printf("%f ", 		C[i*tsize + j*size + k]);
	    }
	    printf("\n");
	}
	printf("\n\n");
    }
}

class Emitter: public adpff::adp_ff_node {
public:
    Emitter(long ntasks):ntasks(ntasks) {}
    void* adp_svc(void*) {
	for(long i=1;i<=ntasks;++i)
	    ff_send_out((void*)i);
	return NULL;
    }
private:
    long ntasks;
};

class Worker: public adpff::adp_ff_node {
public:
    Worker(long size):size(size),tsize(size*size) {
    }
    
    void* adp_svc(void* task) {
	long taskid = (long)(long*)task;
	--taskid;

	double* _A = &A[taskid*tsize];
	double* _B = &B[taskid*tsize];
	double* _C = &C[taskid*tsize];

	for (register int i = 0; i < size; ++i)
	    for (register int j = 0; j < size; ++j) {
		double _Ctmp=0.0;
		for (register int k = 0; k < size; ++k)
		    _Ctmp += _A[(i * size) + k] * _B[(k * size) + j];
		_C[(i * size) + j] = _Ctmp;
	    }
	return GO_ON;
    }
public:
    long size;
    long tsize;
};

int main(int argc, char* argv[]) {

    if (argc<4) {
	printf("use: %s matrix-size streamlen nworkers\n", argv[0]);
	return -1;
    }

    int size      = atoi(argv[1]);
    int streamlen = atoi(argv[2]);
    int nworkers  = atoi(argv[3]);

    int N     = streamlen;
    int tsize = size*size;

    A = (double*)calloc(N, tsize*sizeof(double));
    B = (double*)calloc(N, tsize*sizeof(double));
    C = (double*)calloc(N, tsize*sizeof(double));

    for(int i=0;i<N;++i)
	for(int j=0;j<size;++j)
	    for(int k=0;k<size;++k) {
		A[i*tsize + j*size + k] = (i+j+k)/3.14;
		B[i*tsize + j*size + k] = (i+j*k) + 3.14;
	    }

    Obs obs;
    adpff::AdaptivityParameters ap("demo-fastflow.xml");
    ap.observer = &obs;
    adpff::adp_ff_farm<> farm(ap);
    Emitter E(N);
    farm.add_emitter(&E);
    std::vector<ff_node *> w;
    for(int i=0;i<nworkers;++i) w.push_back(new Worker(size));
    farm.add_workers(w);    

    farm.run_and_wait_end();
    printf("Total time %0.2f (ms)\n", farm.ffTime());
    for(int i=0;i<nworkers;++i)
	printf("Worker %d mean task time %.2f\n", i,
	       diffmsec( ((Worker*)w[i])->getwstoptime(), ((Worker*)w[i])->getwstartime())/(streamlen/nworkers));
    //printf("-------------------\n");Print(N,size);
    return 0;
}
