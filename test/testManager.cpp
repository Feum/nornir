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

std::string getApplicationPath(const std::string& archName){
    return "./" + archName + "/applications/";
}

std::string getBenchPath(const std::string& archName, const std::string& benchName){
    return getApplicationPath(archName) + benchName;
}

std::string getTestPath(const std::string& archName, const std::string& benchName, const std::string& testName){
    return getApplicationPath(archName) + benchName + "/" + testName + "/";
}

std::vector<std::string> getBenchmarks(const std::string& archName){
    return mammut::utils::getFilesNamesInDir(getApplicationPath(archName), false, true);
}

std::vector<std::string> getTestCases(const std::string& archName, const std::string& benchName){
    return mammut::utils::getFilesNamesInDir(getBenchPath(archName, benchName), false, true);
}

typedef struct ApplicationInfo{
    size_t numWorkers;
    size_t numServiceNodes;
}ApplicationInfo;

ApplicationInfo getApplicationInfo(const std::string& archName, const std::string& benchName, const std::string& testName){
    std::vector<std::string> lines = mammut::utils::readFile(getTestPath(archName, benchName, testName) + "info.csv");
    ApplicationInfo info;
    for(auto line : lines){
        if(line.size() && line[0] != '#'){
            auto fields = mammut::utils::split(line, '\t');
            info.numWorkers = mammut::utils::stringToUint(fields[0]);
            info.numServiceNodes = mammut::utils::stringToUint(fields[1]);
        }
    }
    return info;
}

inline std::vector<std::string> getFieldsToCompare(){
    return {"[VirtualCores]", "Workers", "Frequency"};
}

size_t getPosition(const std::string& fieldName, const std::string& headerLine){
    std::vector<std::string> headerFields = mammut::utils::split(headerLine, '\t');
    auto it = std::find(headerFields.begin(), headerFields.end(), fieldName);
    if (it == headerFields.end()){
        throw std::runtime_error("Field " + fieldName + " not present.");
    }else{
        return std::distance(headerFields.begin(), it);
    }
}

TEST(ManagerTest, GlobalTest) {
    // Stored samples were taken with 4 custom fileds.
    EXPECT_EQ(KNARR_MAX_CUSTOM_FIELDS, 4);
    for(std::string arch : getTestingArchitectures()){ // For all architectures
        for(std::string bench : getBenchmarks(arch)){ // For all benchmarks
            for(std::string testcase : getTestCases(arch, bench)){ // And for all testcases on that benchmark
                std::cout << "Running test: " << arch << " - " << bench <<
                             " - " << testcase << std::endl;
                Parameters p = getParameters(arch, getTestPath(arch, bench, testcase) + "parameters.xml");
                // Sampling can be done faster since it is just simulated.
                p.samplingIntervalCalibration = 1;
                p.samplingIntervalSteady = 1;

                // Simulate the execution
                ApplicationInfo info = getApplicationInfo(arch, bench, testcase);
                ManagerTest m(p, info.numWorkers, info.numServiceNodes);
                m.setSimulationParameters(getTestPath(arch, bench, testcase) + "samples.csv");
                m.start();
                m.join();

                /*********************************************************/
                /* Now compare output of the test with the expected one. */
                /*********************************************************/
                std::vector<std::string> expectedLines = mammut::utils::readFile(getTestPath(arch, bench, testcase) + "stats.csv");
                std::vector<std::string> realLines = mammut::utils::readFile("stats.csv");
                // Remove empty lines at end of files
                while(expectedLines.back().compare("") == 0){expectedLines.pop_back();}
                while(realLines.back().compare("") == 0){realLines.pop_back();}
                EXPECT_EQ(expectedLines.size(), realLines.size());
                // Get Fields positions (may be different in different logs).
                // We try to be independent from the position so we keep a map where
                // the first element is the field name and the second is its column
                // position inside the log.
                std::map<std::string, size_t> expectedFields, realFields;
                for(std::string field : getFieldsToCompare()){
                    expectedFields.insert(std::pair<std::string, size_t>(field, getPosition(field, expectedLines[0])));
                    realFields.insert(std::pair<std::string, size_t>(field, getPosition(field, realLines[0])));
                }
                // Now compare the output line by line.
                for(size_t i = 1; i < expectedLines.size(); i++){
                    for(std::string field : getFieldsToCompare()){
                        // Get the fields by name (to be independent from the column position).
                        std::vector<std::string> currentExpectedFields = mammut::utils::split(expectedLines[i], '\t');
                        std::vector<std::string> currentRealFields = mammut::utils::split(realLines[i], '\t');
                        std::string currentExpectedField = currentExpectedFields[expectedFields[field]];
                        std::string currentRealField = currentRealFields[realFields[field]];
                        // For the cores string we need to remove the positions related to service nodes
                        // (the test may be executed with a farm but the validation is done with managertest that
                        // doesn't have service nodes).
                        if(field.compare("[VirtualCores]") == 0){
                            // Remove square brackets
                            currentExpectedField = currentExpectedField.substr(1, currentExpectedField.size() - 2);
                            currentRealField = currentRealField.substr(1, currentRealField.size() - 2);
                            std::vector<std::string> expectedIdentifiers = mammut::utils::split(currentExpectedField, ',');
                            std::vector<std::string> realIdentifiers = mammut::utils::split(currentRealField, ',');
                            EXPECT_EQ(expectedIdentifiers.size() - info.numServiceNodes, realIdentifiers.size());
                            // Compare all the identifiers except those of service nodes.
                            for(size_t i = 0; i < expectedIdentifiers.size() - info.numServiceNodes; i++){
                                EXPECT_EQ(expectedIdentifiers[i], realIdentifiers[i]);
                            }
                        }else{
                            // Compare a generic field.
                            EXPECT_EQ(currentExpectedField, currentRealField);
                        }
                    }
                }
            }
        }
    }
}
