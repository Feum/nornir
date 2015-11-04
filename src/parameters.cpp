/*
 * parameters.hpp
 *
 * Created on: 23/03/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of AdaptiveFastFlow.
 *
 *  AdaptiveFastFlow is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  AdaptiveFastFlow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with AdaptiveFastFlow.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include "parameters.hpp"

namespace adpff{

using namespace std;
using namespace mammut::cpufreq;
using namespace mammut::topology;
using namespace mammut::utils;
using mammut::Communicator;
using mammut::Mammut;
using mammut::utils::enumStrings;

XmlTree::XmlTree(const string& fileName, const string& rootName){
    rapidxml::xml_document<> xmlContent;
    ifstream file(fileName.c_str());
    if(!file.is_open()){
        throw runtime_error("Impossible to read xml file " + fileName);
    }

    string fileContent;
    file.seekg(0, ios::end);
    fileContent.reserve(file.tellg());
    file.seekg(0, ios::beg);
    fileContent.assign((istreambuf_iterator<char>(file)),
                        istreambuf_iterator<char>());
    _fileContentChars = new char[fileContent.size() + 1];
    copy(fileContent.begin(), fileContent.end(), _fileContentChars);
    _fileContentChars[fileContent.size()] = '\0';
    xmlContent.parse<0>(_fileContentChars);
    _root = xmlContent.first_node(rootName.c_str());
}

XmlTree::~XmlTree(){
    delete[] _fileContentChars;
}

void XmlTree::getBool(const char* valueName, bool& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = string(node->value()).compare("true")?false:true;
    }
}

void XmlTree::getInt(const char* valueName, int& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = stringToInt(node->value());
    }
}

void XmlTree::getUint(const char* valueName, uint& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = stringToUint(node->value());
    }
}

void XmlTree::getUlong(const char* valueName, ulong& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = stringToUlong(node->value());
    }
}

void XmlTree::getDouble(const char* valueName, double& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = stringToDouble(node->value());
    }
}

void XmlTree::getString(const char* valueName, string& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        value = node->value();
    }
}

template<typename T>
void XmlTree::getEnum(const char* valueName, T& value){
    rapidxml::xml_node<> *node = NULL;
    node = _root->first_node(valueName);
    if(node){
        stringstream line(node->value());
        line >> enumFromStringInternal(value);
    }
}

void ArchData::loadXml(const string& archFileName){
    XmlTree xc(archFileName, "archData");
    SETVALUE(xc, Double, ticksPerNs);
    SETVALUE(xc, Double, monitoringCost);
    SETVALUE(xc, String, voltageTableFile);
}

void Parameters::setDefault(){
    contractType = CONTRACT_NONE;
    knobWorkers = KNOB_WORKERS_YES;
    knobFrequencies = KNOB_FREQUENCY_YES;
    knobMapping = KNOB_MAPPING_AUTO;
    knobMappingEmitter = KNOB_SNODE_MAPPING_AUTO;
    knobMappingCollector = KNOB_SNODE_MAPPING_AUTO;
    knobHyperthreading = KNOB_HT_AUTO;
    strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_NONE;
    strategyInactiveVirtualCores = STRATEGY_UNUSED_VC_NONE;
    strategyPrediction = STRATEGY_PREDICTION_REGRESSION_LINEAR;
    strategyPredictionErrorPrimary = STRATEGY_PREDICTION_ERROR_CONSTANT;
    strategyPredictionErrorSecondary = STRATEGY_PREDICTION_ERROR_CONSTANT;
    strategySmoothing = STRATEGY_SMOOTHING_EXPONENTIAL;
    strategyCalibration = STRATEGY_CALIBRATION_SOBOL;
    strategyPolling = STRATEGY_POLLING_SLEEP_LATENCY;
    strategyPersistence = STRATEGY_PERSISTENCE_SAMPLES;
    turboBoost = false;
    fastReconfiguration = false;
    migrateCollector = false;
    smoothingFactor = 0;
    persistenceValue = 0;
    samplingInterval = 0;
    underloadThresholdFarm = 80.0;
    overloadThresholdFarm = 90.0;
    underloadThresholdWorker = 80.0;
    overloadThresholdWorker = 90.0;
    requiredBandwidth = 0;
    requiredCompletionTime = 0;
    expectedTasksNumber = 0;
    powerBudget = 0;
    maxPrimaryPredictionError = 5.0;
    maxSecondaryPredictionError = 5.0;
    maxMonitoringOverhead = 1.0;
    observer = NULL;
}

/**
 * Sets the default values for parameters that depends
 * from others.
 */
void Parameters::setDefaultPost(){
    if(!samplingInterval){
        //TODO Se questo sampling interval è molto minore della latenza
        // media di un task potrei settare il sampling interval alla
        // latenza media di un task.
        double msMonitoringCost = (archData.monitoringCost/
                                   archData.ticksPerNs*
                                   0.000001);
        samplingInterval = ceil(msMonitoringCost*
                                    (100.0 - maxMonitoringOverhead));
    }

    if(!smoothingFactor){
        switch(strategySmoothing){
            case STRATEGY_SMOOTHING_MOVING_AVERAGE:{
                smoothingFactor = 10;
            }break;
            case STRATEGY_SMOOTHING_EXPONENTIAL:{
                smoothingFactor = 0.5;
            }break;
        }
    }

    if(!persistenceValue){
        switch(strategyPersistence){
            case STRATEGY_PERSISTENCE_SAMPLES:{
                persistenceValue = 10;
            }break;
            case STRATEGY_PERSISTENCE_TASKS:{
                persistenceValue = 1000;
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

    if(knobFrequencies == KNOB_FREQUENCY_YES &&
       !(isGovernorAvailable(GOVERNOR_USERSPACE) &&
         availableFrequencies.size())){
        knobFrequencies = KNOB_FREQUENCY_NO;
    }

    if(knobFrequencies == KNOB_FREQUENCY_NO){
        for(size_t i = 0; i < virtualCores.size(); i++){
            if(!virtualCores.at(i)->hasFlag("constant_tsc")){
                return VALIDATION_NO_CONSTANT_TSC;
            }
        }
    }else if(knobMapping == KNOB_MAPPING_NO){
        return VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING;
        if(archData.voltageTableFile.empty() ||
           !existsFile(archData.voltageTableFile)){
            return VALIDATION_VOLTAGE_FILE_NEEDED;
        }
    }

    if((fastReconfiguration && !isHighestFrequencySettable()) ||
        knobFrequencies == KNOB_FREQUENCY_NO){
        return VALIDATION_NO_FAST_RECONF;
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobMapping(){
    if(knobMapping == KNOB_MAPPING_AUTO){
        knobMapping = KNOB_MAPPING_LINEAR;
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobSnodeMapping(){
    if(knobMappingEmitter == KNOB_SNODE_MAPPING_AUTO){
        knobMappingEmitter = KNOB_SNODE_MAPPING_ALONE;
    }

    if(knobMappingCollector == KNOB_SNODE_MAPPING_AUTO){
        knobMappingCollector = KNOB_SNODE_MAPPING_ALONE;
    }
    return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobHt(){
    if(knobHyperthreading == KNOB_HT_AUTO){
        knobHyperthreading = KNOB_HT_NO;
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
               strategyPrediction == STRATEGY_PREDICTION_SIMPLE){
                return VALIDATION_WRONG_CONTRACT_PARAMETERS;
            }
        }break;
        default:{
            ;
        }break;
    }

    if(maxPrimaryPredictionError < 0      ||
       maxPrimaryPredictionError > 100.0  ||
       maxSecondaryPredictionError < 0    ||
       maxSecondaryPredictionError > 100.0){
        return VALIDATION_WRONG_CONTRACT_PARAMETERS;
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

template<> char const* enumStrings<KnobConfWorkers>::data[] = {
    "NO",
    "YES"
};

template<> char const* enumStrings<KnobConfFrequencies>::data[] = {
    "NO",
    "YES"
};

template<> char const* enumStrings<KnobConfMapping>::data[] = {
    "NO",
    "AUTO",
    "LINEAR",
    "CACHE_EFFICIENT"
};

template<> char const* enumStrings<KnobConfHyperthreading>::data[] = {
    "NO",
    "AUTO",
    "SOONER",
    "LATER"
};

template<> char const* enumStrings<StrategyUnusedVirtualCores>::data[] = {
    "AUTO",
    "NONE",
    "LOWEST_FREQUENCY",
    "OFF"
};

template<> char const* enumStrings<StrategyPrediction>::data[] = {
    "SIMPLE",
    "REGRESSION_LINEAR"
};

template<> char const* enumStrings<StrategyPredictionError>::data[] = {
    "CONSTANT",
    "COEFFVAR"
};

template<> char const* enumStrings<KnobConfSNodeMapping>::data[] = {
    "NO",
    "AUTO",
    "ALONE",
    "COLLAPSED",
};

template<> char const* enumStrings<StrategySmoothing>::data[] = {
    "MOVING_AVERAGE",
    "EXPONENTIAL"
};

template<> char const* enumStrings<StrategySmoothingFactor>::data[] = {
    "CONST",
    "DYNAMIC"
};

template<> char const* enumStrings<StrategyCalibration>::data[] = {
    "RANDOM",
    "NIEDERREITER",
    "SOBOL",
    "HALTON",
    "HALTON_REVERSE"
};

template<> char const* enumStrings<StrategyPolling>::data[] = {
    "SPINNING",
    "PAUSE",
    "SLEEP_SMALL",
    "SLEEP_LATENCY"
};

template<> char const* enumStrings<StrategyPersistence>::data[] = {
    "SAMPLES",
    "TASKS",
    "VARIATION"
};

void Parameters::loadXml(const string& paramFileName){
    XmlTree xt(paramFileName, "adaptivityParameters");

    SETVALUE(xt, Enum, contractType);
    SETVALUE(xt, Enum, knobMapping);
    SETVALUE(xt, Enum, knobHyperthreading);
    SETVALUE(xt, Enum, knobFrequencies);
    SETVALUE(xt, Enum, strategyUnusedVirtualCores);
    SETVALUE(xt, Enum, strategyInactiveVirtualCores);
    SETVALUE(xt, Enum, strategyPrediction);
    SETVALUE(xt, Enum, strategyPredictionErrorPrimary);
    SETVALUE(xt, Enum, strategyPredictionErrorSecondary);
    SETVALUE(xt, Enum, strategySmoothing);
    SETVALUE(xt, Enum, strategyCalibration);
    SETVALUE(xt, Enum, strategyPolling);
    SETVALUE(xt, Enum, strategyPersistence);
    SETVALUE(xt, Enum, knobMappingEmitter);
    SETVALUE(xt, Enum, knobMappingCollector);

    SETVALUE(xt, Bool, turboBoost);
    SETVALUE(xt, Bool, fastReconfiguration);
    SETVALUE(xt, Double, smoothingFactor);
    SETVALUE(xt, Double, persistenceValue);
    SETVALUE(xt, Uint, samplingInterval);
    SETVALUE(xt, Double, underloadThresholdFarm);
    SETVALUE(xt, Double, overloadThresholdFarm);
    SETVALUE(xt, Double, underloadThresholdWorker);
    SETVALUE(xt, Double, overloadThresholdWorker);
    SETVALUE(xt, Bool, migrateCollector);
    SETVALUE(xt, Double, requiredBandwidth);
    SETVALUE(xt, Uint, requiredCompletionTime);
    SETVALUE(xt, Ulong, expectedTasksNumber);
    SETVALUE(xt, Double, powerBudget);
    SETVALUE(xt, Double, maxPrimaryPredictionError);
    SETVALUE(xt, Double, maxSecondaryPredictionError);
    SETVALUE(xt, Double, maxMonitoringOverhead);
}

Parameters::Parameters(Communicator* const communicator):
      mammut(communicator){
    setDefault();
}

Parameters::Parameters(const string& paramFileName,
                     const string& archFileName,
                     Communicator* const communicator):
      mammut(communicator){
    setDefault();
    archData.loadXml(archFileName);
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

    /** Validate mapping knob. **/
    r = validateKnobMapping();
    if(r != VALIDATION_OK){return r;}

    /** Validate service nodes mapping knob. **/
    r = validateKnobSnodeMapping();
    if(r != VALIDATION_OK){return r;}

    /** Validate hyperthreading knob. **/
    r = validateKnobHt();
    if(r != VALIDATION_OK){return r;}

    /** Validate unused cores strategy. **/
    r = validateUnusedVc(strategyInactiveVirtualCores);
    if(r != VALIDATION_OK){return r;}
    r = validateUnusedVc(strategyUnusedVirtualCores);
    if(r != VALIDATION_OK){return r;}

    /** Validate contract parameters. **/
    r = validateContract();
    if(r != VALIDATION_OK){return r;}

    /** Validate unsupported strategies. **/
    if(knobHyperthreading == KNOB_HT_SOONER ||
       knobMapping == KNOB_MAPPING_CACHE_EFFICIENT){
        throw runtime_error("Not yet supported.");
    }

    return VALIDATION_OK;
}

}

