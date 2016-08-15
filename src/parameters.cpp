/*
 * parameters.hpp
 *
 * Created on: 23/03/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include "parameters.hpp"

#include <cstring>

namespace nornir{

using namespace std;
using namespace mammut::cpufreq;
using namespace mammut::topology;
using namespace mammut::utils;
using mammut::Communicator;
using mammut::Mammut;
using mammut::utils::enumStrings;

void XmlTree::init(const std::string& content, const std::string& rootName){
    rapidxml::xml_document<> xmlContent;
    _fileContentChars = new char[content.size() + 1];
    copy(content.begin(), content.end(), _fileContentChars);
    _fileContentChars[content.size()] = '\0';
    xmlContent.parse<0>(_fileContentChars);
    _root = xmlContent.first_node(rootName.c_str());
}

XmlTree::XmlTree(const string& fileName, const string& rootName){
    ifstream file(fileName.c_str());
    if(!file){
        throw runtime_error("Impossible to read xml file " + fileName);
    }

    string fileContent;
    file.seekg(0, ios::end);
    fileContent.reserve(file.tellg());
    file.seekg(0, ios::beg);
    fileContent.assign((istreambuf_iterator<char>(file)),
                        istreambuf_iterator<char>());
    init(fileContent, rootName);
}

/*XmlTree::XmlTree(const std::string& content, const std::string& rootName){
    init(content, rootName);
}*/

XmlTree::~XmlTree(){
    delete[] _fileContentChars;
}

rapidxml::xml_node<>* XmlTree::getNode(const char* valueName) const{
    rapidxml::xml_node<> *node = _root;
    string vName(valueName);
    size_t i = 0;
    vector<string> tokens = mammut::utils::split(vName, '.');
    if(!tokens.size()){
        return NULL;
    }
    while(i < tokens.size() && node != NULL){
        node = node->first_node(tokens.at(i).c_str());
        i++;
    }
    return node;
}

void XmlTree::getBool(const char* valueName, bool& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = string(node->value()).compare("true")?false:true;
    }
}

void XmlTree::getInt(const char* valueName, int& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = stringToInt(node->value());
    }
}

void XmlTree::getUint(const char* valueName, uint& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = stringToUint(node->value());
    }
}

void XmlTree::getUlong(const char* valueName, ulong& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = stringToUlong(node->value());
    }
}

void XmlTree::getDouble(const char* valueName, double& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = stringToDouble(node->value());
    }
}

void XmlTree::getString(const char* valueName, string& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        value = node->value();
    }
}

void XmlTree::getArrayUint(const char* valueName, std::vector<uint>& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    value.clear();
    if(node){
        vector<string> strValues = mammut::utils::split(node->value(), ':');
        for(size_t i = 0; i < strValues.size(); i++){
            value.push_back(stringToUint(strValues.at(i)));
        }
    }
}


template<typename T>
void XmlTree::getEnum(const char* valueName, T& value){
    rapidxml::xml_node<> *node = getNode(valueName);
    if(node){
        stringstream line(node->value());
        line >> enumFromStringInternal(value);
    }
}

void ArchData::loadXml(const string& archFileName){
    XmlTree xc(archFileName, "archData");
    SETVALUE(xc, Double, ticksPerNs);
    SETVALUE(xc, Double, monitoringCost);
}

void Parameters::setDefault(){
    contractType = CONTRACT_NONE;
    triggerQBlocking = TRIGGER_Q_BLOCKING_NO;
    strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_SAME;
    strategySelection = STRATEGY_SELECTION_LEARNING;
    strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_USL;
    strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR;
    strategyExploration = STRATEGY_EXPLORATION_HALTON;
    strategySmoothing = STRATEGY_SMOOTHING_EXPONENTIAL;
    strategyPolling = STRATEGY_POLLING_SLEEP_SMALL;
    strategyPersistence = STRATEGY_PERSISTENCE_SAMPLES;
    strategyCoresChange = STRATEGY_CORES_RETHREADING;
    knobCoresEnabled = true;
    knobMappingEnabled = true;
    knobFrequencyEnabled = true;
    turboBoost = false;
    fastReconfiguration = true;
    migrateCollector = false;
    smoothingFactor = 0;
    persistenceValue = 0;
    samplingIntervalCalibration = 100;
    samplingIntervalSteady = 1000;
    steadyThreshold = 4;
    minTasksPerSample = 0;
    underloadThresholdFarm = 80.0;
    overloadThresholdFarm = 90.0;
    underloadThresholdWorker = 80.0;
    overloadThresholdWorker = 90.0;
    requiredBandwidth = 0;
    requiredCompletionTime = 0;
    expectedTasksNumber = 0;
    synchronousWorkers = false;
    powerBudget = 0;
    maxCalibrationTime = 0;
    maxPerformancePredictionError = 10.0;
    maxPowerPredictionError = 5.0;
    regressionAging = 0;
    maxMonitoringOverhead = 1.0;
    thresholdQBlocking = -1;
    tolerableSamples = 0;
    qSize = 1;
    conservativeValue = 0;
    isolateManager = false;
    statsReconfiguration = false;

    mishra.applicationName = "";
    mishra.namesData = "";
    mishra.bandwidthData = "";
    mishra.powerData = "";

    dataflow.orderedProcessing = false;
    dataflow.orderedOutput = false;
    dataflow.maxGraphs = 1000;
    dataflow.maxInterpreters = 0;

    observer = NULL;
}

uint Parameters::getLowOverheadSamplingInterval() const{
    //TODO Se questo sampling interval è molto minore della latenza
    // media di un task potrei settare il sampling interval alla
    // latenza media di un task.

    // We are setting the sampling interval such that the latency
    // due to the communication between manager and node is the
    // maxMonitoringOverhead% of the interval.
    double msMonitoringCost = (archData.monitoringCost/
                               archData.ticksPerNs*
                               0.000001);
    return ceil(msMonitoringCost*(100.0 - maxMonitoringOverhead));
}

/**
 * Sets the default values for parameters that depends
 * from others.
 */
void Parameters::setDefaultPost(){
    if(!samplingIntervalCalibration){
        samplingIntervalCalibration = getLowOverheadSamplingInterval();
    }

    if(!samplingIntervalSteady){
        samplingIntervalSteady = getLowOverheadSamplingInterval();
    }

    if(!smoothingFactor){
        switch(strategySmoothing){
            case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
                smoothingFactor = 10;
            }break;
            case STRATEGY_SMOOTHING_EXPONENTIAL:{
                smoothingFactor = 0.1;
            }break;
        }
    }

    if(!persistenceValue){
        switch(strategyPersistence){
            case STRATEGY_PERSISTENCE_SAMPLES:{
                persistenceValue = 1;
            }break;
            case STRATEGY_PERSISTENCE_VARIATION:{
                persistenceValue = 5;
            }break;
        }
    }
}

bool Parameters::isGovernorAvailable(Governor g){
    return mammut.getInstanceCpuFreq()->isGovernorAvailable(g);
}

vector<Frequency> Parameters::getAvailableFrequencies(){
    vector<Frequency> frequencies;
    CpuFreq* cpuFreq = mammut.getInstanceCpuFreq();
    if(cpuFreq){
        vector<Domain*> fDomains = cpuFreq->getDomains();
        if(fDomains.size()){
            frequencies = fDomains.front()->getAvailableFrequencies();
        }
    }
    return frequencies;
}

bool Parameters::isUnusedVcOffAvailable(){
    vector<VirtualCore*> vc = mammut.getInstanceTopology()->
                              getVirtualCores();
    for(size_t i = 0; i < vc.size(); i++){
        if(vc.at(i)->isHotPluggable()){
            return true;
        }
    }
    return false;
}

bool Parameters::isFrequencySettable(){
    vector<Frequency> frequencies = getAvailableFrequencies();
    return  isGovernorAvailable(GOVERNOR_USERSPACE) && frequencies.size();

}

bool Parameters::isLowestFrequencySettable(){
    return isGovernorAvailable(GOVERNOR_POWERSAVE) ||
           isFrequencySettable();
}

bool Parameters::isHighestFrequencySettable(){
    return isGovernorAvailable(GOVERNOR_PERFORMANCE) ||
           isFrequencySettable();
}

ParametersValidation Parameters::validateUnusedVc(StrategyUnusedVirtualCores& s){
    switch(s){
        case STRATEGY_UNUSED_VC_OFF:{
            if(!isUnusedVcOffAvailable()){
                return VALIDATION_UNUSED_VC_NO_OFF;
            }
        }break;
        case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
            if(!isLowestFrequencySettable()){
                return VALIDATION_UNUSED_VC_NO_FREQUENCIES;
            }
        }break;
        case STRATEGY_UNUSED_VC_AUTO:{
            if(isUnusedVcOffAvailable()){
                s = STRATEGY_UNUSED_VC_OFF;
            }else if(isLowestFrequencySettable()){
                s = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY;
            }else{
                s = STRATEGY_UNUSED_VC_SAME;
            }
        }break;
        default:
            break;
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobFrequencies(){
    vector<Frequency> availableFrequencies = getAvailableFrequencies();
    vector<VirtualCore*> virtualCores;
    virtualCores = mammut.getInstanceTopology()->getVirtualCores();

    if(knobFrequencyEnabled && !(isGovernorAvailable(GOVERNOR_USERSPACE) &&
         availableFrequencies.size())){
        knobFrequencyEnabled = false;
    }

    for(size_t i = 0; i < virtualCores.size(); i++){
        if(!virtualCores.at(i)->hasFlag("constant_tsc")){
            return VALIDATION_NO_CONSTANT_TSC;
        }
    }

    if(knobFrequencyEnabled && !knobMappingEnabled){
        return VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING;
    }

    if(fastReconfiguration &&
       (!isHighestFrequencySettable() ||
         strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_NONE)){
        fastReconfiguration = false;
    }

    if(mammut.getInstanceCpuFreq()->isBoostingSupported()){
        if(turboBoost){
            mammut.getInstanceCpuFreq()->enableBoosting();
        }else{
            mammut.getInstanceCpuFreq()->disableBoosting();
        }
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateTriggers(){
    if(triggerQBlocking == TRIGGER_Q_BLOCKING_AUTO &&
       thresholdQBlocking == -1){
        return VALIDATION_NO_BLOCKING_THRESHOLD;
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateContract(){
    switch(contractType){
        case CONTRACT_PERF_UTILIZATION:{
            if(underloadThresholdFarm > overloadThresholdFarm     ||
               underloadThresholdWorker > overloadThresholdWorker ||
               underloadThresholdFarm < 0                         ||
               overloadThresholdFarm > 100                        ||
               underloadThresholdWorker < 0                       ||
               overloadThresholdWorker > 100){
                return VALIDATION_WRONG_CONTRACT_PARAMETERS;
            }
        }break;
        case CONTRACT_PERF_BANDWIDTH:{
            if(requiredBandwidth <= 0){
                return VALIDATION_WRONG_CONTRACT_PARAMETERS;
            }
        }break;
        case CONTRACT_PERF_COMPLETION_TIME:{
            if(!expectedTasksNumber || !requiredCompletionTime){
                return VALIDATION_WRONG_CONTRACT_PARAMETERS;
            }
        }break;
        case CONTRACT_POWER_BUDGET:{
            if(powerBudget <= 0 ||
               strategySelection == STRATEGY_SELECTION_ANALYTICAL ||
               strategySelection == STRATEGY_SELECTION_LIMARTINEZ){
                return VALIDATION_WRONG_CONTRACT_PARAMETERS;
            }
        }break;
        default:{
            ;
        }break;
    }

    if(maxCalibrationTime == 0 &&
       (maxPerformancePredictionError <= 0      ||
        maxPerformancePredictionError > 100.0   ||
        maxPowerPredictionError <= 0    ||
        maxPowerPredictionError > 100.0)){
        return VALIDATION_WRONG_CONTRACT_PARAMETERS;
    }

    return VALIDATION_OK;
}

ParametersValidation Parameters::validateSelector(){
    if(strategySelection == STRATEGY_SELECTION_MISHRA &&
       (mishra.bandwidthData.compare("") == 0 ||
        mishra.powerData.compare("") == 0 ||
        mishra.applicationName.compare("") == 0 ||
        mishra.namesData.compare("") == 0)){
        return VALIDATION_NO_MISHRA_PARAMETERS;
    }
    return VALIDATION_OK;
}


ParametersValidation Parameters::validatePredictor(){
    if(strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_MISHRA &&
         (mishra.bandwidthData.compare("") == 0 ||
          mishra.powerData.compare("") == 0 ||
          mishra.applicationName.compare("") == 0 ||
          mishra.namesData.compare("") == 0)){
            return VALIDATION_NO_MISHRA_PARAMETERS;
    }
    return VALIDATION_OK;
}

template<> char const* enumStrings<ContractType>::data[] = {
    "NONE",
    "PERF_UTILIZATION",
    "PERF_BANDWIDTH",
    "PERF_COMPLETION_TIME",
    "POWER_BUDGET"
};

template<> char const* enumStrings<TriggerConfQBlocking>::data[] = {
    "NO",
    "YES",
    "AUTO"
};

template<> char const* enumStrings<StrategyUnusedVirtualCores>::data[] = {
    "AUTO",
    "SAME",
    "NONE",
    "LOWEST_FREQUENCY",
    "OFF"
};

template<> char const* enumStrings<StrategySelection>::data[] = {
    "LEARNING",
    "ANALYTICAL",
    "FULLSEARCH",
    "LIMARTINEZ",
    "MISHRA"
};

template<> char const* enumStrings<StrategyPredictionPerformance>::data[] = {
    "AMDAHL",
    "AMDAHL_MAPPING",
    "USL",
    "USL_MAPPING",
    "MISHRA"
};

template<> char const* enumStrings<StrategyPredictionPower>::data[] = {
    "LINEAR",
    "LINEAR_MAPPING",
    "MISHRA"
};

template<> char const* enumStrings<StrategyExploration>::data[] = {
    "RANDOM",
    "NIEDERREITER",
    "SOBOL",
    "HALTON",
    "HALTON_REVERSE"
};

template<> char const* enumStrings<StrategySmoothing>::data[] = {
    "MOVING_AVERAGE",
    "EXPONENTIAL"
};

template<> char const* enumStrings<StrategySmoothingFactor>::data[] = {
    "CONST",
    "DYNAMIC"
};

template<> char const* enumStrings<StrategyPolling>::data[] = {
    "SPINNING",
    "PAUSE",
    "SLEEP_SMALL",
    "SLEEP_LATENCY"
};

template<> char const* enumStrings<StrategyPersistence>::data[] = {
    "SAMPLES",
    "VARIATION"
};

template<> char const* enumStrings<StrategyCoresChange>::data[] = {
    "RETHREADING",
    "REMAPPING"
};

void Parameters::loadXml(const string& paramFileName){
    XmlTree xt(paramFileName, "adaptivityParameters");

    SETVALUE(xt, Enum, contractType);
    SETVALUE(xt, Enum, strategyUnusedVirtualCores);
    SETVALUE(xt, Enum, strategySelection);
    SETVALUE(xt, Enum, strategyPredictionPerformance);
    SETVALUE(xt, Enum, strategyPredictionPower);
    SETVALUE(xt, Enum, strategyExploration);
    SETVALUE(xt, Enum, strategySmoothing);
    SETVALUE(xt, Enum, strategyPolling);
    SETVALUE(xt, Enum, strategyPersistence);
    SETVALUE(xt, Enum, strategyCoresChange);
    SETVALUE(xt, Enum, triggerQBlocking);
    SETVALUE(xt, Bool, knobCoresEnabled);
    SETVALUE(xt, Bool, knobMappingEnabled);
    SETVALUE(xt, Bool, knobFrequencyEnabled);

    SETVALUE(xt, Bool, turboBoost);
    SETVALUE(xt, Bool, fastReconfiguration);
    SETVALUE(xt, Double, smoothingFactor);
    SETVALUE(xt, Double, persistenceValue);
    SETVALUE(xt, Uint, samplingIntervalCalibration);
    SETVALUE(xt, Uint, samplingIntervalSteady);
    SETVALUE(xt, Uint, steadyThreshold);
    SETVALUE(xt, Uint, minTasksPerSample);
    SETVALUE(xt, Double, underloadThresholdFarm);
    SETVALUE(xt, Double, overloadThresholdFarm);
    SETVALUE(xt, Double, underloadThresholdWorker);
    SETVALUE(xt, Double, overloadThresholdWorker);
    SETVALUE(xt, Bool, migrateCollector);
    SETVALUE(xt, Double, requiredBandwidth);
    SETVALUE(xt, Uint, requiredCompletionTime);
    SETVALUE(xt, Ulong, expectedTasksNumber);
    SETVALUE(xt, Bool, synchronousWorkers);
    SETVALUE(xt, Double, powerBudget);
    SETVALUE(xt, Double, maxCalibrationTime);
    SETVALUE(xt, Double, maxPerformancePredictionError);
    SETVALUE(xt, Double, maxPowerPredictionError);
    SETVALUE(xt, Uint, regressionAging);
    SETVALUE(xt, Double, maxMonitoringOverhead);
    SETVALUE(xt, Double, thresholdQBlocking);
    SETVALUE(xt, Uint, tolerableSamples);
    SETVALUE(xt, Ulong, qSize);
    SETVALUE(xt, Double, conservativeValue);
    SETVALUE(xt, ArrayUint, disallowedNumCores);
    SETVALUE(xt, Bool, isolateManager);
    SETVALUE(xt, Bool, statsReconfiguration);

    SETVALUE(xt, String, mishra.applicationName);
    SETVALUE(xt, String, mishra.namesData);
    SETVALUE(xt, String, mishra.bandwidthData);
    SETVALUE(xt, String, mishra.powerData);

    SETVALUE(xt, Bool, dataflow.orderedProcessing);
    SETVALUE(xt, Bool, dataflow.orderedOutput);
    SETVALUE(xt, Uint, dataflow.maxGraphs);
    SETVALUE(xt, Uint, dataflow.maxInterpreters);
}

Parameters::Parameters(Communicator* const communicator):
      mammut(communicator){
    setDefault();
}

#define CONFIGURATION_VERSION "\"1.0.0\""
#define CONFPATH_LEN_MAX 512
#define CONFFILE_VERSION "/nornir/version.csv"
#define CONFFILE_ARCH "/nornir/archdata.xml"
#define CONFFILE_VOLTAGE "/nornir/voltage.csv"

Parameters::Parameters(const string& paramFileName,
                       Communicator* const communicator):
      mammut(communicator){
    setDefault();

    /** Retrieving archdata.xml configuration file. **/
    char* confHome_c = getenv("XDG_CONFIG_DIRS");
    vector<string> confHomes;
    if(!confHome_c || strcmp(confHome_c, "") == 0){
        confHomes.push_back(string("/etc/xdg"));
    }else{
        confHomes = split(string(confHome_c), ':');
    }

    size_t i = 0;
    bool found = false;
    while(i < confHomes.size() && !found){
        string confFileArch = confHomes.at(i) + string(CONFFILE_ARCH);
        string confFileVoltage = confHomes.at(i) + string(CONFFILE_VOLTAGE);
        string confFileVersion = confHomes.at(i) + string(CONFFILE_VERSION);

        if(existsFile(confFileArch) &&
           existsFile(confFileVoltage) &&
           existsFile(confFileVersion) &&
           readFirstLineFromFile(confFileVersion).compare(CONFIGURATION_VERSION) == 0){
            archData.loadXml(confFileArch);
            loadVoltageTable(archData.voltageTable, confFileVoltage);
            found = true;
        }
        i++;
    }

    if(!found){
        throw runtime_error("Impossible to find configuration files. Please run 'sudo make microbench' from the nornir root folder.");
    }

    /** Loading parameters. **/
    loadXml(paramFileName);
}

Parameters::~Parameters(){
    ;
}

ParametersValidation Parameters::validate(){
    ParametersValidation r = VALIDATION_OK;

    setDefaultPost();

    /** Validate frequency knob. **/
    r = validateKnobFrequencies();
    if(r != VALIDATION_OK){return r;}

    /** Validate triggers. **/
    r = validateTriggers();
    if(r != VALIDATION_OK){return r;}

    /** Validate unused cores strategy. **/
    r = validateUnusedVc(strategyUnusedVirtualCores);
    if(r != VALIDATION_OK){return r;}

    /** Validate contract parameters. **/
    r = validateContract();
    if(r != VALIDATION_OK){return r;}

    /** Validate selectors. **/
    r = validateSelector();
    if(r != VALIDATION_OK){return r;}

    /** Validate predictors. **/
    r = validatePredictor();
    if(r != VALIDATION_OK){return r;}

    return VALIDATION_OK;
}

}

