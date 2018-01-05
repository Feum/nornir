/**
 *  Different tests on parallel for.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "../src/nornir.hpp"
#include "gtest/gtest.h"

using namespace nornir;


void runTest(int startloop, int endloop, int step, int chunksize){
    int nworkers = 4;
    std::vector<uint> v;
    nornir::Parameters p = getParameters("repara");
    nornir::ParallelFor pf(nworkers, &p);
    for(size_t i = 0; i < 10; i++){
        v.resize(endloop - startloop, 0);

        pf.parallel_for(startloop, endloop, step, chunksize, 
        [&v](long long int idx, long long int id){
            v[idx] = idx;
        });

        for(int j = startloop; j < endloop; j += step){
            if((j - startloop) % step == 0){
                EXPECT_EQ(v[j], (uint) j);
            }else{
                EXPECT_EQ(v[j], (uint) 0);
            }
        }
        v.clear();
    }    
}

TEST(ParallelForTest, Simple){
    runTest(0, 100, 1, 0);
}

TEST(ParallelForTest, Step3){
    runTest(0, 100, 3, 0);
}

TEST(ParallelForTest, DynamicChunk1){
    runTest(0, 100, 3, 1);
}

TEST(ParallelForTest, DynamicChunk3){
    runTest(0, 100, 3, 3);
}