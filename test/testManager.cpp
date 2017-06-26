/**
 *  Different tests on Manager.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h> 
#include <linux/limits.h>
#include <time.h>
#include "../src/manager.hpp"
#include "gtest/gtest.h"

using namespace nornir;

typedef struct ApplicationInfo{
    size_t numServiceNodes;
}ApplicationInfo;

std::map<std:string, ApplicationInfo> applicationInfos;

TEST(SelectorLearnerTest, SimpleTest) {
    #if 0
    ApplicationInfo ai;
    
    ai.numServiceNodes = 2;
    applicationInfos["pbzip2"] = ai;
    //TODO Parametric execution wrt benchmarkname, architecture and infos.
    #endif

    // Stored samples are taken with 4 custom fileds.
    EXPECT_EQ(KNARR_MAX_CUSTOM_FIELDS, 4);
    std::vector<std::string> testsFiles;
    testsFiles.push_back("./repara/applications/pbzip2/watts_30_learning_usl_");

    for(std::string test : testsFiles){
        std::cout << "Running test: " << test << std::endl;
        Parameters p = getParametersRepara(test + "parameters.xml");
        // Sampling can be done faster since it is just simulated.
        //p.samplingIntervalCalibration = 1;
        //p.samplingIntervalSteady = 1;
        
        p.knobHyperthreadingEnabled = false;
        nornir::Observer o;
        p.observer = &o;
        ManagerTest m(p, 2);
        m.setSimulationParameters(test + "samples.csv");
        m.start();
        m.join();
    }
}
