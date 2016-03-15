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
import sys


fields = {}

fields['Percentile'] = 0
fields['PrimaryRequired'] = 1
fields['PrimaryLossCnt'] = 2
fields['PrimaryAvg'] = 3
fields['PrimaryStddev'] = 4
fields['PrimaryLossAvg'] = 5
fields['PrimaryLossStddev'] = 6
fields['SecondaryOptimal'] = 7
fields['SecondaryAvg'] = 8
fields['SecondaryStddev'] = 9
fields['SecondaryLossAvg'] = 10
fields['SecondaryLossStddev'] = 11
fields['CalibrationStepsAvg'] = 12
fields['CalibrationStepsStddev'] = 13
fields['CalibrationTimeAvg'] = 14
fields['CalibrationTimeStddev'] = 15
fields['CalibrationTimePercAvg'] = 16
fields['CalibrationTimePercStddev'] = 17
fields['CalibrationTaksAvg'] = 18
fields['CalibrationTasksStddev'] = 19
fields['CalibrationTasksPercAvg'] = 20
fields['CalibrationTasksPercStddev'] = 21
fields['ReconfigurationTimeWorkersAvg'] = 22
fields['ReconfigurationTimeWorkersStddev'] = 23
fields['ReconfigurationTimeFrequencyAvg'] = 24
fields['ReconfigurationTimeFrequencyStddev'] = 25
fields['ReconfigurationTimeTotalAvg'] = 26
fields['ReconfigurationTimeTotalStddev'] = 27

alternatives = ('REGRESSION_LINEAR_RANDOM', 'REGRESSION_LINEAR_HALTON', 'REGRESSION_LINEAR_HALTON_FAST')
alternativesReconfTime = ('REGRESSION_LINEAR_HALTON', 'REGRESSION_LINEAR_HALTON_FAST')
alternativesMandelPerf = ('LIMARTINEZ', 'REGRESSION_LINEAR_HALTON_FAST', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE10', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE20', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE30')
alternativesMandelPower = ('LIMARTINEZ', 'REGRESSION_LINEAR_HALTON_FAST', 'REGRESSION_LINEAR_HALTON_FAST_AGING3', 'REGRESSION_LINEAR_HALTON_FAST_AGING6', 'REGRESSION_LINEAR_HALTON_FAST_AGING10')



benchmarks = ('blackscholes', 'canneal', 'pbzip2', 'simple_mandelbrot', 'videoprocessing')
#contracts = ('PERF_COMPLETION_TIME', 'POWER_BUDGET')
contracts = ('POWER_BUDGET')

iterations = 5

def getBound(parameters):
    bound = -1
    parametersFile = open(parameters)
    for line in parametersFile:
        if 'powerBudget' in line:
            bound = float(line.split('<')[1].split('>')[1])
    parametersFile.close()
    return bound

for bench in benchmarks:
    print ""
    if 'pbzip2' in bench:
        iterations = 1
    else:
        iterations = 5
    for c in xrange(0,1):
        for alt in alternatives:
            averageTimePerc = []
            averageWattsPerc = []
            for boundperc in xrange(10, 100, 20):
                iterationsTimePerc = []
                iterationsWattsPerc = []
                iterationsJoulesOverhead = []
                for i in xrange(0, iterations):
                    previousTs = 0
                    totalJoules = 0
                    prefix = 'RESULTS/' + bench + '/POWER_BUDGET_' + alt + '/' + str(boundperc) + '.' + str(i)
                    #print prefix
                    summary = open(prefix + '.summary.csv')
                    time = 0
                    realWatts = 0
                    realJoules = 0
                    for line in summary:
                        time = line.split('\t')[2]
                        realWatts = line.split('\t')[0]
                    realJoules = float(time) * float(realWatts)
                    summary.close()

                    stats = open(prefix + '.stats.csv')
                    bound = getBound(prefix + '.parameters.xml')
                    violationSecs = 0
                    violationJoules = 0

                    buckets = [0] * int(math.ceil(float(time) + 1))
                    firstline = 1
                    for line in stats:
                        if firstline:
                            firstline = 0
                            continue
                        currentTs = float(line.split('\t')[0]) / 1000.0
                        currentWatts = float(line.split('\t')[10])
                        currentJoules = currentWatts * (currentTs - previousTs)
                        totalJoules = totalJoules + currentJoules

                        previousBucket = int(math.floor(previousTs))
                        currentBucket = int(math.floor(currentTs))
                        if previousBucket != currentBucket:
                            if currentBucket > previousBucket + 1:
                                intbuckets = (currentBucket - previousBucket + 1) - 2
                                for midb in xrange(previousBucket + 1, currentBucket):
                                    buckets[midb] += currentWatts
                                currentJoules -= intbuckets*currentWatts

                            buckets[previousBucket] += (math.ceil(previousTs) - previousTs)*currentWatts
                            buckets[currentBucket] += (currentTs - math.floor(currentTs))*currentWatts
                            #buckets[previousBucket] += (currentBucket - previousTs)*currentJoules
                            #buckets[currentBucket] += (currentTs - currentBucket)*currentJoules
                        else:
                            buckets[currentBucket] += currentJoules

                        previousTs = currentTs

                    secondsOut = 0
                    joulesOut = 0
                    totalSec = 0
                    wattsOut = 0 
                    for b in buckets:
                        if b:
                            totalSec += 1
                        if b > bound:
                            secondsOut += 1
                            joulesOut += b
                    if secondsOut:
                        wattsOut = joulesOut / secondsOut

                    iterationsTimePerc.append((float(secondsOut)/totalSec) * 100.0)
                    if wattsOut:
                        iterationsWattsPerc.append(float((wattsOut - bound) / bound) * 100.0)
                    else:
                        iterationsWattsPerc.append(0)
                    iterationsJoulesOverhead.append(realJoules - totalJoules)

                    #print "TotalSecs: " + str(totalSec)
                    #print "SecsOut: " + str(secondsOut)
                    #print "JoulesOut: " + str(joulesOut)
                    #print "PercTimeOut: " + str((float(secondsOut)/totalSec) * 100.0)
                    #print "PercWattsOut: " + str(((wattsOut - bound) / bound) * 100.0)
                #print str(np.average(iterationsJoulesOverhead))
                averageTimePerc.append(np.average(iterationsTimePerc))
                averageWattsPerc.append(np.average(iterationsWattsPerc))
            print bench + " " + alt + " AvgTime: " + str(np.average(averageTimePerc))
            print bench + " " + alt + " AvgWatts: " + str(np.average(averageWattsPerc))


