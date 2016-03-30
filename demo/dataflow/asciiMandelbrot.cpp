#include <iostream>
#include <stdlib.h>
#include "../../src/dataflow/manager.hpp"

namespace mandel{

using namespace nornir::dataflow;

class MandelInputStream:public nornir::dataflow::InputStream{
private:
	int x,y;
	bool eos,testAsync;
	clock_t t;
public:
	inline MandelInputStream(bool asy):x(-39),y(-39),eos(false),testAsync(asy),t(clock()){;}

       inline nornir::dataflow::StreamElem* next(){
		if(testAsync && clock()-t<=2000)
			return NULL;
		else{
			if(testAsync) t=clock();
			if(x==39&&y==38)
				eos=true;
			if(x==39){
				x=-39;
				++y;
				ArrayWrapper<float> *t=new ArrayWrapper<float>(2);
				t->set(0,-100000);
				return t;
			}
			ArrayWrapper<float> *t=new ArrayWrapper<float>(2);
			t->set(0,x/40.0f);
			t->set(1,y/40.0f);
			++x;
			return t;
		}
	}

	inline bool hasNext(){return !eos;}

};

class MandelOutputStream:public nornir::dataflow::OutputStream{
public:
        void put(nornir::dataflow::StreamElem* a){
		Wrapper<float> *dt=(Wrapper<float>*) a;
		float x=dt->get();
		if(x==0)
			std::cout << "*";
		else if(x==-1)
			std::cout << std::endl;
		else
			std::cout << " ";
		delete a;
	}
};

int BAILOUT = 16;
int MAX_ITERATIONS;



/**Function taken from http://www.timestretch.com/FractalBenchmark.html#67b4f5a3200c7b7e900c38ff21321741 **/
Wrapper<float>* iterateTask(ArrayWrapper<float>* t){
	float f1=t->get(1);
	float f2=t->get(0);
	delete t;
	if(f2==-100000)
		return new Wrapper<float>(-1);

	float cr = f1-0.5f;
	float ci = f2;
	float zi = 0.0f;
	float zr = 0.0f;
	int i = 0;
		while (true) {
		i++;
		float temp = zr * zi;
		float zr2 = zr * zr;
		float zi2 = zi * zi;
		zr = zr2 - zi2 + cr;
		zi = temp + temp + ci;
		if (zi2 + zr2 > BAILOUT)
			return new Wrapper<float>(i);

		if (i > MAX_ITERATIONS)
			return new Wrapper<float>(0);

	}
}

void exec(int mit,int nWorkers,int groupSize, bool async){
	MAX_ITERATIONS=mit;
	time_t start=time(NULL);
	MandelInputStream inp(async);
	MandelOutputStream out;
	nornir::dataflow::Farm* f = nornir::dataflow::createStandardFarm<ArrayWrapper<float>,Wrapper<float>,iterateTask>();
	nornir::dataflow::Manager m(f,&inp,&out,nWorkers,groupSize,false);
	m.exec();
	std::cout << "Esecuzione con "<<nWorkers<<" workers in: " << time(NULL) - start << std::endl;
	delete f;
}


}

int main(int argc,char** argv){
	using namespace mandel;
	if(argc<5){
		std::cerr << "Usage:" << std::endl;
		std::cerr << argv[0] <<" maxIterations parDegree groupSize testAsync" << std::endl;
		return -1;
	}else{
		int nElem=atoi(argv[1]);
		int parDegree=atoi(argv[2]);
		int groupSize=atoi(argv[3]);
		bool async=atoi(argv[4]);
		exec(nElem,parDegree,groupSize,async);
		return 0;
	}
}

