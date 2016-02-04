import math
import subprocess
import shlex
import socket
import os
import os.path
import shutil
import argparse 
from subprocess import Popen, PIPE, STDOUT
import numpy as np

RESULTS_FILE = 'REPARA_results.csv'
powersList = []
timesList = []

iterations = 10 

def getLastLine(fileName):
    fh = open(fileName, "r")
    for line in fh:
        pass
    last = line
    fh.close()
    return last


def run(contractType, fieldName, fieldValue, percentile, itnum):
    parametersFile = open("parameters.xml", "w")
    parametersFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    parametersFile.write("<adaptivityParameters>\n")

    parametersFile.write("<qSize>1</qSize>\n")
    parametersFile.write("<statsReconfiguration>true</statsReconfiguration>\n")
    parametersFile.write("<samplingIntervalCalibration>500</samplingIntervalCalibration>\n")
    parametersFile.write("<samplingIntervalSteady>2000</samplingIntervalSteady>\n")
    parametersFile.write("<maxPrimaryPredictionError>10.0</maxPrimaryPredictionError>\n")
    parametersFile.write("<smoothingFactor>0.01</smoothingFactor>\n")
    parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
    parametersFile.write("<minTasksPerSample>1</minTasksPerSample>\n")
    parametersFile.write("<persistenceValue>1</persistenceValue>\n")
    parametersFile.write("<strategyPolling>SLEEP_SMALL</strategyPolling>\n")
    parametersFile.write("<strategyPrediction>" + args.prediction + "</strategyPrediction>\n")
    parametersFile.write("<contractType>" + contractType + "</contractType>\n")
    parametersFile.write("<" + fieldName + ">" + str(fieldValue) + "</" + fieldName + ">\n")
    parametersFile.write("</adaptivityParameters>\n")
    parametersFile.close()

    run = "./pbzip2_ff -f -k -p22 /home/desensi/enwiki-20151201-abstract.xml"
    process = subprocess.Popen(shlex.split(run), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = process.communicate()

    logfile = open("log.csv", "w")
    logfile.write(out)
    logfile.write(err)
    logfile.close()

    time = -1
    watts = -1
    calibrationSteps = -1
    calibrationTime = -1
    calibrationTimePerc = -1
    calibrationTasks = -1
    calibrationTasksPerc = -1
    reconfigurationTimeAvg = -1
    reconfigurationTimeStddev = -1
    try:
        result = getLastLine("summary.csv")
        time = float(result.split('\t')[2])
        watts = float(result.split('\t')[0])

        calibrationSteps = float(result.split('\t')[3])
        calibrationTime = float(result.split('\t')[4])
        calibrationTimePerc = float(result.split('\t')[5])
        calibrationTasks = float(result.split('\t')[6])
        calibrationTasksPerc = float(result.split('\t')[7])
        reconfigurationTimeAvg = float(result.split('\t')[8])
        reconfigurationTimeStddev = float(result.split('\t')[9])
    except:
        print "nothing"

    dir_results = "nodir"
    if contractType == 'PERF_COMPLETION_TIME':
        dir_results = perfdir
    elif contractType == 'POWER_BUDGET':
        dir_results = powerdir

    for outfile in ["calibration.csv", "stats.csv", "summary.csv", "parameters.xml", "log.csv"]:
        if os.path.isfile(outfile):
            os.rename(outfile, dir_results + "/" + str(percentile) + "." + str(itnum) + "." + outfile)

    # Returns data
    return time, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeAvg

def loadPerfPowerData():
    fh = open(RESULTS_FILE, "r")
    for line in fh:
        #Workers Frequency Time Watts
        if line[0] != '#':
            fields = line.split("\t")
            timesList.append(float(fields[2]))
            powersList.append(float(fields[3]))

    fh.close()

def getOptimalTimeBound(time):
    fh = open(RESULTS_FILE, "r")
    bestPower = 99999999999
    for line in fh:
        #Workers Frequency Time Watts
        if line[0] != '#':
            fields = line.split("\t")
            curtime = float(fields[2])
            curpower = float(fields[3])
            if curpower < bestPower and curtime <= float(time):
                bestPower = curpower
        
    fh.close()
    return bestPower

def getOptimalPowerBound(power):
    fh = open(RESULTS_FILE, "r")
    bestTime = 9999999999
    for line in fh:
        #Workers Frequency Time Watts
        if line[0] != '#':
            fields = line.split("\t")
            curtime = float(fields[2])
            curpower = float(fields[3])
            if curtime < bestTime and curpower <= float(power):
                bestTime = curtime
    fh.close()
    return bestTime

############################################################################################################

parser = argparse.ArgumentParser(description='Runs the application to check the accuracy of the reconfiguration.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p', '--prediction', help='Prediction strategy.', required=True)
parser.add_argument('-c', '--contract', help='Contract type.', required=True)
args = parser.parse_args()

loadPerfPowerData()

if args.contract == "PERF_COMPLETION_TIME":
    perfdir = 'PERF_COMPLETION_TIME_' + args.prediction
    shutil.rmtree(perfdir, ignore_errors=True)
    os.makedirs(perfdir)
    perfFile = open(perfdir + "/results.csv", "w")
    outFile = perfFile
elif args.contract == "POWER_BUDGET":
    powerdir = 'POWER_BUDGET_' + args.prediction
    shutil.rmtree(powerdir, ignore_errors=True)
    os.makedirs(powerdir)
    powerFile = open(powerdir + "/results.csv", "w")
    outFile = powerFile

outFile.write("#Percentile\tPrimaryRequired\tPrimaryAvg\tPrimaryStddev\tPrimaryLossAvg\tPrimaryLossStddev\t" \
              "SecondaryOptimal\tSecondaryAvg\tSecondaryStddev\tSecondaryLossAvg\tSecondaryLossStddev\t" \
              "CalibrationStepsAvg\tCalibrationStepsStddev\t" \
              "CalibrationTimeAvg\tCalibrationTimeStddev\tCalibrationTimePercAvg\tCalibrationTimePercStddev\t"
              "CalibrationTaksAvg\tCalibrationTasksStddev\tCalibrationTasksPercAvg\tCalibrationTasksPercStddev\t"
              "ReconfigurationTimeAvg\tReconfigurationTimeStddev\n")

totalSecondaryLosses = []
totalCalibrationSteps = []
totalCalibrationTime = []
totalCalibrationTimePerc = []
totalCalibrationTasks = []
totalCalibrationTasksPerc = []
totalReconfigurationTime = []

for p in xrange(10, 100, 20):
    cts = []
    wattses = []
    opts = []
    calibrationsSteps = []
    calibrationsTime = []
    calibrationsTimePerc = []
    calibrationsTasks = []
    calibrationsTasksPerc = []
    primaryLosses = []
    secondaryLosses = []
    reconfigurationTimes = []
    avgPrimary = 0
    stddevPrimary = 0
    avgSecondary = 0
    stddevSecondary = 0
    avgLoss = 0
    stddevLoss = 0
    avgCalibrationSteps = 0
    stddevCalibrationSteps = 0
    avgCalibrationTime = 0
    stddevCalibrationTime = 0
    avgCalibrationTimePerc = 0
    stddevCalibrationTimePerc = 0
    avgCalibrationTasks = 0
    stddevCalibrationTasks = 0
    avgCalibrationTasksPerc = 0
    stddevCalibrationTasksPerc = 0
    avgReconfigurationTime = 0
    stddevReconfigurationTime = 0
    primaryReq = 0
    target = 0

    for i in xrange(0, iterations):
        ct = -1
        if args.contract == "PERF_COMPLETION_TIME":
            ############################
            # targetTime = np.percentile(np.array(timesList), p)
            ############################
            minTime = min(timesList)
            maxTime = max(timesList)
            targetTime = minTime + (maxTime - minTime)*(p/100.0)  
            target = targetTime
            while ct == -1:
                ct, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeAvg = run("PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, p, i)
            opt = getOptimalTimeBound(targetTime)
            primaryReq = targetTime
            primaryLoss = ((ct - targetTime) / targetTime ) * 100.0
            secondaryLoss = ((watts - opt) / opt) * 100.0
        elif args.contract == "POWER_BUDGET":
            #############################
            # targetPower = np.percentile(np.array(powersList), p)
            #############################
            minPower = min(powersList)
            maxPower = max(powersList)
            targetPower = minPower + (maxPower - minPower)*(p/100.0)
            target = targetPower
            while ct == -1:
                ct, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeAvg = run("POWER_BUDGET", "powerBudget", targetPower, p, i)
            opt = getOptimalPowerBound(targetPower)
            primaryReq = targetPower
            primaryLoss = ((watts - targetPower) / targetPower) * 100.0
            secondaryLoss = ((ct - opt) / opt) * 100.0
        cts.append(ct)
        wattses.append(watts)
        opts.append(opt)
        calibrationsSteps.append(calibrationSteps)
        calibrationsTime.append(calibrationTime)
        calibrationsTimePerc.append(calibrationTimePerc)
        calibrationsTasks.append(calibrationTasks)
        calibrationsTasksPerc.append(calibrationTasksPerc)
        primaryLosses.append(primaryLoss)
        secondaryLosses.append(secondaryLoss)
        reconfigurationTimes.append(reconfigurationTimeAvg)

    if args.contract == "PERF_COMPLETION_TIME":
        avgPrimary = np.average(cts)
        stddevPrimary = np.std(cts)
        avgSecondary = np.average(wattses)
        stddevSecondary = np.std(wattses)
    elif args.contract == "POWER_BUDGET":    
        avgPrimary = np.average(wattses)
        stddevPrimary = np.std(wattses)
        avgSecondary = np.average(cts)
        stddevSecondary = np.std(cts)
    
    avgPrimaryLoss = np.average(primaryLosses)
    stddevPrimaryLoss = np.std(primaryLosses)
    avgSecondaryLoss = np.average(secondaryLosses)
    stddevSecondaryLoss = np.std(secondaryLosses)

    avgCalibrationSteps = np.average(calibrationsSteps)
    stddevCalibrationSteps = np.std(calibrationsSteps)
    avgCalibrationTime = np.average(calibrationsTime)
    stddevCalibrationTime = np.std(calibrationsTime)
    avgCalibrationTimePerc = np.average(calibrationsTimePerc)
    stddevCalibrationTimePerc = np.std(calibrationsTimePerc)
    avgCalibrationTasks = np.average(calibrationsTasks)
    stddevCalibrationTasks = np.std(calibrationsTasks)
    avgCalibrationTasksPerc = np.average(calibrationsTasksPerc)
    stddevCalibrationTasksPerc = np.std(calibrationsTasksPerc)
    avgReconfigurationTime = np.average(reconfigurationTimes)
    stddevReconfigurationTime = np.std(reconfigurationTimes)

    totalSecondaryLosses.append(avgSecondaryLoss)
    totalCalibrationSteps.append(avgCalibrationSteps)
    totalCalibrationTime.append(avgCalibrationTime)
    totalCalibrationTimePerc.append(avgCalibrationTimePerc)
    totalCalibrationTasks.append(avgCalibrationTasks)
    totalCalibrationTasksPerc.append(avgCalibrationTasksPerc)
    totalReconfigurationTime.append(avgReconfigurationTime)

    outFile.write(str(p) + "\t" + str(target) + "\t" + str(avgPrimary) + "\t" + str(stddevPrimary) + "\t" + str(avgPrimaryLoss) + "\t" + str(stddevPrimaryLoss) + "\t" + 
                                  str(opt) + "\t" + str(avgSecondary) + "\t" + str(stddevSecondary) + "\t" + str(avgSecondaryLoss) + "\t" + str(stddevSecondaryLoss) + "\t" + 
                                  str(avgCalibrationSteps) + "\t" + str(stddevCalibrationSteps) + "\t" + 
                                  str(avgCalibrationTime) + "\t" + str(stddevCalibrationTime) + "\t" + str(avgCalibrationTimePerc) + "\t" + str(stddevCalibrationTimePerc) + "\t" + 
                                  str(avgCalibrationTasks) + "\t" + str(stddevCalibrationTasks) + "\t" + str(avgCalibrationTasksPerc) + "\t" + str(stddevCalibrationTasksPerc) + "\t" +
                                  str(avgReconfigurationTime) + "\t" + str(stddevReconfigurationTime) + "\n")
    outFile.flush()
    os.fsync(outFile.fileno())


### Now print the average over all the possible targets
outFile.write("=================================================================================================================================\n")
outFile.write("=                                                         TOTAL RESULTS                                                         =\n")
outFile.write("=================================================================================================================================\n")
outFile.write("AvgSecondaryLoss\tStddevSecondaryLoss\t" \
              "AvgCalibrationSteps\tStddevCalibrationSteps\t" \
              "AvgCalibrationTime\tStddevCalibrationTime\tAvgCalibrationTimePerc\tStddevCalibrationTimePerc\t" \
              "AvgCalibrationTasks\tStddevCalibrationTasks\tAvgCalibrationTasksPerc\tStddevCalibrationTasksPerc\t" \
              "AvgReconfigurationTime\tStddevReconfigurationTime\n")

outFile.write(str(np.average(totalSecondaryLosses)) + "\t" + str(np.std(totalSecondaryLosses)) + "\t" +
              str(np.average(totalCalibrationSteps))+ "\t" + str(np.std(totalCalibrationSteps)) + "\t" +
              str(np.average(totalCalibrationTime))+ "\t" + str(np.std(totalCalibrationTime)) + "\t" + str(np.average(totalCalibrationTimePerc))+ "\t" + str(np.std(totalCalibrationTimePerc)) + "\t" +
              str(np.average(totalCalibrationTasks))+ "\t" + str(np.std(totalCalibrationTasks)) + "\t" + str(np.average(totalCalibrationTasksPerc))+ "\t" + str(np.std(totalCalibrationTasksPerc)) + "\t" +
              str(np.average(totalReconfigurationTime))+ "\t" + str(np.std(totalReconfigurationTime)) + "\n")

if args.contract == "PERF_COMPLETION_TIME":
    perfFile.close()
elif args.contract == "POWER_BUDGET":
    powerFile.close()
