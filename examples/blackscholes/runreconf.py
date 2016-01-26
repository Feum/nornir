import math
import subprocess
import shlex
import socket
import os
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


def run(contractType, fieldName, fieldValue, percentile):
    parametersFile = open("parameters.xml", "w")
    parametersFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    parametersFile.write("<adaptivityParameters>\n")

    parametersFile.write("<qSize>4</qSize>\n")
    parametersFile.write("<samplingIntervalCalibration>100</samplingIntervalCalibration>\n")
    parametersFile.write("<samplingIntervalSteady>1000</samplingIntervalSteady>\n")
    parametersFile.write("<smoothingFactor>0.1</smoothingFactor>\n")
    parametersFile.write("<strategyPersistence>TASKS</strategyPersistence>\n")
    parametersFile.write("<persistenceValue>1000</persistenceValue>\n")
    parametersFile.write("<strategyPolling>SLEEP_LATENCY</strategyPolling>\n")
    parametersFile.write("<strategyPrediction>" + args.prediction + "</strategyPrediction>\n")
    parametersFile.write("<contractType>" + contractType + "</contractType>\n")
    parametersFile.write("<" + fieldName + ">" + str(fieldValue) + "</" + fieldName + ">\n")
    parametersFile.write("</adaptivityParameters>\n")
    parametersFile.close()

    run = "./blackscholes 23 inputs/in_10M.txt tmp.txt"
    process = subprocess.Popen(shlex.split(run), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = process.communicate()

    time = -1
    watts = -1
    calibration = -1
    try:
        result = getLastLine("summary.csv")
        time = float(result.split('\t')[2])
        watts = float(result.split('\t')[0])
        calibration = float(result.split('\t')[3])
    except:
        print "nothing"

    dir_results = "nodir"
    if contractType == 'PERF_COMPLETION_TIME':
        dir_results = perfdir
    elif contractType == 'POWER_BUDGET':
        dir_results = powerdir

    for outfile in ["calibration.csv", "stats.csv", "summary.csv", "parameters.xml"]:
        os.rename(outfile, dir_results + "/" + str(percentile) + "." + outfile)

    # Returns completion time and watts
    return time, watts, calibration

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

outFile.write("#Fraction\tPrimaryRequired\tPrimaryAvg\tPrimaryStddev\tSecondaryOptimal\tSecondaryAvg\tSecondaryStddev\tLossAvg\tLossStddev\tCalibrationAvg\tCalibrationStddev\n")

for p in xrange(10, 110, 10):
    cts = []
    wattses = []
    opts = []
    calibrations = []
    avgPrimary = 0
    stddevPrimary = 0
    avgSecondary = 0
    stddevSecondary = 0
    avgLoss = 0
    stddevLoss = 0
    avgCalibration = 0
    stddevCalibration = 0
    primaryReq = 0

    for i in xrange(0, iterations):
        if args.contract == "PERF_COMPLETION_TIME":
            targetTime = np.percentile(np.array(timesList), p)
            ct, watts, calibration = run("PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, p)
            opt = getOptimalTimeBound(targetTime)
            primaryReq = targetTime
            loss = ((watts - opt) / opt) * 100.0
        elif args.contract == "POWER_BUDGET":
            targetPower = np.percentile(np.array(powersList), p)
            ct, watts, calibration = run("POWER_BUDGET", "powerBudget", targetPower, p)
            opt = getOptimalPowerBound(targetPower)
            primaryReq = targetPower
            loss = ((ct - opt) / opt) * 100.0
        cts.append(ct)
        wattses.append(watts)
        opts.append(opt)
        calibrations.append(calibration)

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

    avgLoss = np.average(losses)
    stddevLoss = np.std(losses)
    avgCalibration = np.average(calibrations)
    stddevCalibration = np.std(calibrations)

    outFile.write(str(p) + "\t" + str(avgPrimary) + "\t" + str(stddevPrimary) + "\t" + str(opt) + "\t" + str(avgSecondary) + "\t" + str(stddevSecondary) + "\t" + str(avgLoss) + "\t" + str(stddevLoss) + "\t" + str(avgCalibration) + "\t" + str(stddevCalibration) + "\n")
    outFile.flush()
    os.fsync(outFile.fileno())

if args.contract == "PERF_COMPLETION_TIME":
    perfFile.close()
elif args.contract == "POWER_BUDGET":
    powerFile.close()
