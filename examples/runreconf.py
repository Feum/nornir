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

iterations = 5

def getLastLine(fileName):
    fh = open(fileName, "r")
    for line in fh:
        pass
    last = line
    fh.close()
    return last

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def run(bench, contractType, fieldName, fieldValue, percentile, dir_results, calstrategy, fastreconf, aging, conservative, knobworkers, itnum):
    with cd(bench):
        parametersFile = open("parameters.xml", "w")
        parametersFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        parametersFile.write("<adaptivityParameters>\n")
        parametersFile.write("<strategyPolling>SLEEP_SMALL</strategyPolling>\n")
        parametersFile.write("<statsReconfiguration>true</statsReconfiguration>\n")
        parametersFile.write("<strategyPrediction>" + args.prediction + "</strategyPrediction>\n")
        parametersFile.write("<contractType>" + contractType + "</contractType>\n")
        parametersFile.write("<" + fieldName + ">" + str(fieldValue) + "</" + fieldName + ">\n")
        
        maxPerformanceError = 10.0
        if bench == 'simple_mandelbrot':
            maxPerformanceError = 20.0

        if contractType == "PERF_COMPLETION_TIME":
            parametersFile.write("<maxPrimaryPredictionError>" + str(maxPerformanceError) + "</maxPrimaryPredictionError>\n")
        else:
            parametersFile.write("<maxSecondaryPredictionError>" + str(maxPerformanceError) + "</maxSecondaryPredictionError>\n")

        if fastreconf:
            print "Using fast reconf."
            parametersFile.write("<fastReconfiguration>true</fastReconfiguration>\n")

        if calstrategy is not None:
            print "Using strategy " + calstrategy
            parametersFile.write("<strategyCalibration>" + calstrategy + "</strategyCalibration>\n")

        if knobworkers is not None:
            parametersFile.write("<knobWorkers>" + knobworkers + "</knobWorkers>\n")

        if bench == 'blackscholes':
            run = "./blackscholes 23 inputs/in_10M.txt tmp.txt"
            parametersFile.write("<qSize>4</qSize>\n")
            parametersFile.write("<samplingIntervalCalibration>100</samplingIntervalCalibration>\n")
            parametersFile.write("<samplingIntervalSteady>1000</samplingIntervalSteady>\n")
            parametersFile.write("<smoothingFactor>0.1</smoothingFactor>\n")
            parametersFile.write("<minTasksPerSample>500</minTasksPerSample>\n")
            parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
            parametersFile.write("<persistenceValue>1</persistenceValue>\n")
        elif bench == 'pbzip2':
            run = "./pbzip2_ff -f -k -p22 /home/desensi/enwiki-20151201-abstract.xml"
            parametersFile.write("<qSize>1</qSize>\n")
            parametersFile.write("<samplingIntervalCalibration>500</samplingIntervalCalibration>\n")
            parametersFile.write("<samplingIntervalSteady>2000</samplingIntervalSteady>\n")
            parametersFile.write("<smoothingFactor>0.001</smoothingFactor>\n")
            parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
            parametersFile.write("<persistenceValue>1</persistenceValue>\n")
        elif bench == 'canneal':
            run = "./canneal 23 15000 2000 2500000.nets 6000"
            parametersFile.write("<qSize>1</qSize>\n")
            parametersFile.write("<samplingIntervalCalibration>10</samplingIntervalCalibration>\n")
            parametersFile.write("<samplingIntervalSteady>1000</samplingIntervalSteady>\n")
            parametersFile.write("<smoothingFactor>0.01</smoothingFactor>\n")
            parametersFile.write("<minTasksPerSample>10</minTasksPerSample>\n")
            parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
            parametersFile.write("<persistenceValue>1</persistenceValue>\n")
        elif bench == 'videoprocessing':
            run = "./video 0 1 0 22 VIRAT/960X540.mp4 VIRAT/480X270.mp4 VIRAT/480X270.mp4 VIRAT/960X540.mp4"
            parametersFile.write("<qSize>1</qSize>\n")
            parametersFile.write("<samplingIntervalCalibration>10</samplingIntervalCalibration>\n")
            parametersFile.write("<samplingIntervalSteady>1000</samplingIntervalSteady>\n")
            parametersFile.write("<minTasksPerSample>5</minTasksPerSample>\n")
            parametersFile.write("<expectedTasksNumber>64790</expectedTasksNumber>\n")
            parametersFile.write("<smoothingFactor>0.7</smoothingFactor>\n")
            parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
            parametersFile.write("<persistenceValue>1</persistenceValue>\n")
        elif bench == 'simple_mandelbrot':
            run = "./mandel_ff 8192 4096 1 22"
            parametersFile.write("<qSize>4</qSize>\n")
            parametersFile.write("<samplingIntervalCalibration>50</samplingIntervalCalibration>\n")
            parametersFile.write("<samplingIntervalSteady>1000</samplingIntervalSteady>\n")
            parametersFile.write("<smoothingFactor>0.9</smoothingFactor>\n")
            if aging is not None:
                parametersFile.write("<regressionAging>" + str(aging) + "</regressionAging>\n")
            parametersFile.write("<minTasksPerSample>4</minTasksPerSample>\n")
            parametersFile.write("<strategyPersistence>SAMPLES</strategyPersistence>\n")
            parametersFile.write("<persistenceValue>1</persistenceValue>\n")
            if conservative is not None:
                parametersFile.write("<conservativeValue>" + str(conservative) + "</conservativeValue>\n")

        parametersFile.write("</adaptivityParameters>\n")
        parametersFile.close()

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
        reconfigurationTimeWorkersAvg = -1
        reconfigurationTimeWorkersStddev = -1
        reconfigurationTimeFrequencyAvg = -1
        reconfigurationTimeFrequencyStddev = -1
        reconfigurationTimeTotalAvg = -1
        reconfigurationTimeTotalStddev = -1
        try:
            result = getLastLine("summary.csv")
            time = float(result.split('\t')[2])
            watts = float(result.split('\t')[0])

            calibrationSteps = float(result.split('\t')[3])
            calibrationTime = float(result.split('\t')[4])
            calibrationTimePerc = float(result.split('\t')[5])
            calibrationTasks = float(result.split('\t')[6])
            calibrationTasksPerc = float(result.split('\t')[7])
            reconfigurationTimeWorkersAvg = float(result.split('\t')[8])
            reconfigurationTimeWorkersStddev = float(result.split('\t')[9])
            reconfigurationTimeFrequencyAvg = float(result.split('\t')[10])
            reconfigurationTimeFrequencyStddev = float(result.split('\t')[11])
            reconfigurationTimeTotalAvg = float(result.split('\t')[14])
            reconfigurationTimeTotalStddev = float(result.split('\t')[15])
        except:
            print "nothing"

        for outfile in ["calibration.csv", "stats.csv", "summary.csv", "parameters.xml", "log.csv"]:
            if os.path.isfile(outfile):
                os.rename(outfile, dir_results + "/" + str(percentile) + "." + str(itnum) + "." + outfile)

    # Returns data
    return time, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeWorkersAvg, reconfigurationTimeFrequencyAvg, reconfigurationTimeTotalAvg

def loadPerfPowerData(bench):
    with cd(bench):
        fh = open(RESULTS_FILE, "r")
        for line in fh:
            #Workers Frequency Time Watts
            if line[0] != '#':
                fields = line.split("\t")
                timesList.append(float(fields[2]))
                powersList.append(float(fields[3]))

        fh.close()

def getOptimalTimeBound(time, bench):
    with cd(bench):
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

def getOptimalPowerBound(power, bench):
    with cd(bench):
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
parser.add_argument('-b', '--benchmark', help='Benchmark.', required=True)
parser.add_argument('-p', '--prediction', help='Prediction strategy.', required=True)
parser.add_argument('-c', '--contract', help='Contract type.', required=True)
parser.add_argument('-k', '--calstrategy', help='Calibration strategy.', required=False, default='HALTON')
parser.add_argument('-f', '--fastreconf', help='Fast reconfiguration.', required=False,  action='store_true')
parser.add_argument('-a', '--aging', help='Aging factor.', required=False)
parser.add_argument('-o', '--conservative', help='Conservative value.', required=False)
parser.add_argument('-w', '--knobworkers', help='Knob Workers.', required=False)

args = parser.parse_args()

loadPerfPowerData(args.benchmark)

rdir = args.contract + '_'
rdir += args.prediction 
if args.prediction == 'REGRESSION_LINEAR':
    rdir += '_' + args.calstrategy
if args.fastreconf:
    rdir += '_FAST'
if args.knobworkers:
    rdir += '_' + args.knobworkers
if args.aging:
    rdir += '_AGING' + str(args.aging)
if args.conservative:
    rdir += '_CONSERVATIVE' + str(args.conservative)

with cd(args.benchmark):
    shutil.rmtree(rdir, ignore_errors=True)
    os.makedirs(rdir)
    outFile = open(rdir + "/results.csv", "w")

outFile.write("#Percentile\tPrimaryRequired\tPrimaryLossCnt\tPrimaryAvg\tPrimaryStddev\tPrimaryLossAvg\tPrimaryLossStddev\t" \
              "SecondaryOptimal\tSecondaryAvg\tSecondaryStddev\tSecondaryLossAvg\tSecondaryLossStddev\t" \
              "CalibrationStepsAvg\tCalibrationStepsStddev\t" \
              "CalibrationTimeAvg\tCalibrationTimeStddev\tCalibrationTimePercAvg\tCalibrationTimePercStddev\t"
              "CalibrationTaksAvg\tCalibrationTasksStddev\tCalibrationTasksPercAvg\tCalibrationTasksPercStddev\t"
              "ReconfigurationTimeWorkersAvg\tReconfigurationTimeWorkersStddev\t"
              "ReconfigurationTimeFrequencyAvg\tReconfigurationTimeFrequencyStddev\t"
              "ReconfigurationTimeTotalAvg\tReconfigurationTimeTotalStddev\n")

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
    reconfigurationTimesWorkers = []
    reconfigurationTimesFrequency = []
    reconfigurationTimesTotal = []
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
    avgReconfigurationTimeWorkers = 0
    stddevReconfigurationTimeWorkers = 0
    avgReconfigurationTimeFrequency = 0
    stddevReconfigurationTimeFrequency = 0
    avgReconfigurationTimeTotal = 0
    stddevReconfigurationTimeTotal = 0
    primaryReq = 0
    target = 0
    
    if args.benchmark == 'pbzip2':
        iterations = 1

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
                ct, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeWorkersAvg, reconfigurationTimeFrequencyAvg, reconfigurationTimeTotalAvg = run(args.benchmark, "PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, p, rdir, args.calstrategy, args.fastreconf, args.aging, args.conservative, args.knobworkers, i)
            opt = getOptimalTimeBound(targetTime, args.benchmark)
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
                ct, watts, calibrationSteps, calibrationTime, calibrationTimePerc, calibrationTasks, calibrationTasksPerc, reconfigurationTimeWorkersAvg, reconfigurationTimeFrequencyAvg, reconfigurationTimeTotalAvg = run(args.benchmark, "POWER_BUDGET", "powerBudget", targetPower, p, rdir, args.calstrategy, args.fastreconf, args.aging, args.conservative, args.knobworkers, i)
            opt = getOptimalPowerBound(targetPower, args.benchmark)
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
        if primaryLoss > 0:
            primaryLosses.append(primaryLoss)
        else:
            secondaryLosses.append(secondaryLoss)
        reconfigurationTimesWorkers.append(reconfigurationTimeWorkersAvg)
        reconfigurationTimesFrequency.append(reconfigurationTimeFrequencyAvg)
        reconfigurationTimesTotal.append(reconfigurationTimeTotalAvg)

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
    
    if len(primaryLosses):
        avgPrimaryLoss = np.average(primaryLosses)
        stddevPrimaryLoss = np.std(primaryLosses)
    else:
        avgPrimaryLoss = 0
        stddevPrimaryLoss = 0

    if len(secondaryLosses):
        avgSecondaryLoss = np.average(secondaryLosses)
        stddevSecondaryLoss = np.std(secondaryLosses)
    else:
        avgSecondaryLoss = 0
        stddevSecondaryLoss = 0

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
    avgReconfigurationTimeWorkers = np.average(reconfigurationTimesWorkers)
    stddevReconfigurationTimeWorkers = np.std(reconfigurationTimesWorkers)
    avgReconfigurationTimeFrequency = np.average(reconfigurationTimesFrequency)
    stddevReconfigurationTimeFrequency = np.std(reconfigurationTimesFrequency)
    avgReconfigurationTimeTotal = np.average(reconfigurationTimesTotal)
    stddevReconfigurationTimeTotal = np.std(reconfigurationTimesTotal)

    outFile.write(str(p) + "\t" + str(target) + "\t" + str(len(primaryLosses)) + "\t" + str(avgPrimary) + "\t" + str(stddevPrimary) + "\t" + str(avgPrimaryLoss) + "\t" + str(stddevPrimaryLoss) + "\t" + 
                                  str(opt) + "\t" + str(avgSecondary) + "\t" + str(stddevSecondary) + "\t" + str(avgSecondaryLoss) + "\t" + str(stddevSecondaryLoss) + "\t" + 
                                  str(avgCalibrationSteps) + "\t" + str(stddevCalibrationSteps) + "\t" + 
                                  str(avgCalibrationTime) + "\t" + str(stddevCalibrationTime) + "\t" + str(avgCalibrationTimePerc) + "\t" + str(stddevCalibrationTimePerc) + "\t" + 
                                  str(avgCalibrationTasks) + "\t" + str(stddevCalibrationTasks) + "\t" + str(avgCalibrationTasksPerc) + "\t" + str(stddevCalibrationTasksPerc) + "\t" +
                                  str(avgReconfigurationTimeWorkers) + "\t" + str(stddevReconfigurationTimeWorkers) + "\t" +
                                  str(avgReconfigurationTimeFrequency) + "\t" + str(stddevReconfigurationTimeFrequency) + "\t" + 
                                  str(avgReconfigurationTimeTotal) + "\t" + str(stddevReconfigurationTimeTotal) + "\n")
    outFile.flush()
    os.fsync(outFile.fileno())

outFile.close()

