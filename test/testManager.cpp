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
using std::vector;
using std::string;

string getApplicationPath(const string& archName){
    return "./" + archName + "/applications/";
}

string getBenchPath(const string& archName, const string& benchName){
    return getApplicationPath(archName) + benchName;
}

string getTestPath(const string& archName,
                   const string& benchName,
                   const string& testName){
    return getApplicationPath(archName) + benchName + "/" + testName + "/";
}

vector<string> getBenchmarks(const string& archName){
    return mammut::utils::getFilesNamesInDir(getApplicationPath(archName), false, true);
}

vector<string> getTestCases(const string& archName, const string& benchName){
    return mammut::utils::getFilesNamesInDir(getBenchPath(archName, benchName), false, true);
}

typedef struct ApplicationInfo{
    size_t numWorkers;
    size_t numServiceNodes;
}ApplicationInfo;

ApplicationInfo getApplicationInfo(const string& archName,
                                   const string& benchName,
                                   const string& testName){
    vector<string> lines = mammut::utils::readFile(getTestPath(archName, benchName, testName) + "info.csv");
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

inline vector<string> getFieldsToCompare(){
    return {"[VirtualCores]", "Workers", "Frequency"};
}

size_t getPosition(const string& fieldName, const string& headerLine){
    vector<string> headerFields = mammut::utils::split(headerLine, '\t');
    auto it = std::find(headerFields.begin(), headerFields.end(), fieldName);
    if (it == headerFields.end()){
        throw std::runtime_error("Field " + fieldName + " not present.");
    }else{
        return std::distance(headerFields.begin(), it);
    }
}

/**
 * We get as input the monitored data, and we verify that the
 * algorithm take the same reconfiguration decision we have in
 * the stored sample files.
 */
TEST(ManagerTest, GlobalTest) {
    // Stored samples were taken with 4 custom fileds.
    EXPECT_EQ(KNARR_MAX_CUSTOM_FIELDS, 4);
    EXPECT_EQ(MAX_RHO, 96);
    for(string arch : getTestingArchitectures()){ // For all architectures
        for(string bench : getBenchmarks(arch)){ // For all benchmarks
            for(string testcase : getTestCases(arch, bench)){ // And for all testcases on that benchmark
                std::cout << "Running test: " << arch << " - " << bench <<
                             " - " << testcase << std::endl;
                Parameters  p = getParameters(arch, getTestPath(arch, bench, testcase) + "parameters.xml");
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
                vector<string> expectedLines = mammut::utils::readFile(getTestPath(arch, bench, testcase) + "stats.csv");
                vector<string> realLines = mammut::utils::readFile("stats.csv");
                // Remove empty lines at end of files
                while(expectedLines.back().compare("") == 0){expectedLines.pop_back();}
                while(realLines.back().compare("") == 0){realLines.pop_back();}
                EXPECT_EQ(expectedLines.size(), realLines.size());
                // Get Fields positions (may be different in different logs).
                // We try to be independent from the position so we keep a map where
                // the first element is the field name and the second is its column
                // position inside the log.
                std::map<string, size_t> expectedFields, realFields;
                for(string field : getFieldsToCompare()){
                    expectedFields.insert(std::pair<string, size_t>(field, getPosition(field, expectedLines[0])));
                    realFields.insert(std::pair<string, size_t>(field, getPosition(field, realLines[0])));
                }
                // Now compare the output line by line.
                for(size_t i = 1; i < expectedLines.size(); i++){
                    for(string field : getFieldsToCompare()){
                        // Get the fields by name (to be independent from the column position).
                        vector<string> currentExpectedFields = mammut::utils::split(expectedLines[i], '\t');
                        vector<string> currentRealFields = mammut::utils::split(realLines[i], '\t');
                        string currentExpectedField = currentExpectedFields[expectedFields[field]];
                        string currentRealField = currentRealFields[realFields[field]];
                        // For the cores string we need to remove the positions related to service nodes
                        // (the test may be executed with a farm but the validation is done with managertest that
                        // doesn't have service nodes).
                        if(field.compare("[VirtualCores]") == 0){
                            // Remove square brackets
                            currentExpectedField = currentExpectedField.substr(1, currentExpectedField.size() - 2);
                            currentRealField = currentRealField.substr(1, currentRealField.size() - 2);
                            vector<string> expectedIdentifiers = mammut::utils::split(currentExpectedField, ',');
                            vector<string> realIdentifiers = mammut::utils::split(currentRealField, ',');
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
