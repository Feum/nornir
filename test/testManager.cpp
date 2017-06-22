/**
 *  Different tests on Manager.
 **/
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "../src/manager.hpp"
#include "gtest/gtest.h"

using namespace nornir;

Parameters getParametersRepara(std::string fileName){
    //p.strategyUnusedVirtualCores // For knobFrequency
    //p.isolateManager //For knobMapping
    //p.knobHyperthreadingEnabled // For knobCores
    //p.disallowedNumCores // For knobCores
    Parameters p(fileName);
    mammut::SimulationParameters simulationParameters;
    simulationParameters.sysfsRootPrefix = "../src/external/mammut/test/archs/repara";
    p.mammut.setSimulationParameters(simulationParameters);
    return p;
}

TEST(SelectorLearnerTest, SimpleTest) {
    // Stored samples are taken with 4 custom fileds.
    EXPECT_EQ(KNARR_MAX_CUSTOM_FIELDS, 4);
    std::vector<std::string> testsFiles;
    testsFiles.push_back("./repara/pbzip2/watts_30_learning_usl_");

    for(std::string test : testsFiles){
        std::cout << "Running test: " << test << std::endl;
        Parameters p = getParametersRepara(test + "parameters.xml");
        ManagerInstrumented m("dummy", p);
        m.setSimulationParameters(test + "samples.csv");
        m.start();
        m.join();
    }
}
