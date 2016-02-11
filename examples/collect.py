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


percFieldPrimaryLoss = 4
percFieldSecondaryLoss = 9

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

alternatives = ('LIMARTINEZ', 'REGRESSION_LINEAR_RANDOM', 'REGRESSION_LINEAR_HALTON', 'REGRESSION_LINEAR_HALTON_FAST')
alternativesReconfTime = ('REGRESSION_LINEAR_HALTON', 'REGRESSION_LINEAR_HALTON_FAST')
alternativesMandelPerf = ('LIMARTINEZ', 'REGRESSION_LINEAR_HALTON_FAST', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE10', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE20', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE30')
alternativesMandelPower = ('LIMARTINEZ', 'REGRESSION_LINEAR_HALTON_FAST', 'REGRESSION_LINEAR_HALTON_FAST_AGING3', 'REGRESSION_LINEAR_HALTON_FAST_AGING6', 'REGRESSION_LINEAR_HALTON_FAST_AGING10')


def printField(bench, contract, alt, field):
    try:
        fh = open(bench + "/" + contract + "_" + alt + "/results.csv")
        fh.readline() # Skip header
        values = []

        for i in xrange(10, 100, 20):
            fieldsLine = fh.readline()
            fieldValue = float(fieldsLine.split('\t')[fields[field]])
            if field != 'PrimaryLossAvg' or fieldValue != 0:
                values.append(fieldValue)

        if field == 'PrimaryLossCnt':
            sys.stdout.write(str(np.sum(values)) + '\t0\t')
        else:
            if len(values):
                avg = np.average(values)
                std = np.std(values)
            else:
                avg = 0
                std = 0
            sys.stdout.write(str(avg) + '\t' + str(std) + '\t')
        fh.close()
    except:
        sys.stdout.write('N.D.\tN.D.\t')

def getAlt(alt):
    alt = alt.replace('REGRESSION_LINEAR_HALTON', 'RLH')
    alt = alt.replace('REGRESSION_LINEAR', 'RL')
    alt = alt.replace('LIMARTINEZ', 'LM')
    return alt

#############################################################################################
benchmarks = ('blackscholes', 'canneal', 'pbzip2', 'simple_mandelbrot', 'videoprocessing')
contracts = ('PERF_COMPLETION_TIME', 'POWER_BUDGET')
fields 

parser = argparse.ArgumentParser(description='Runs the application to check the accuracy of the reconfiguration.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-b', '--benchmark', help='Benchmark.', required=False)
parser.add_argument('-c', '--contract', help='Contract type.', required=True)
parser.add_argument('-f', '--field', help='Field type (PrimaryLossCnt, PrimaryLossAvg, SecondaryLossAvg, CalibrationTimePercAvg, ReconfigurationTimeWorkersAvg).', required=True)
args = parser.parse_args()


if args.benchmark is not None:
    benchmarks = ([args.benchmark])
    if 'simple_mandelbrot' in args.benchmark:
        if 'PERF' in args.contract:
            alternatives = alternativesMandelPerf
        else:
            alternatives = alternativesMandelPower

if args.contract is not None:
    contracts = ([args.contract])

sys.stdout.write('#Bench\t')

if 'ReconfigurationTime' in args.field:
    for alt in alternativesReconfTime:
        sys.stdout.write(getAlt(alt) + '_AVG\t' + getAlt(alt) + '_STDDEV\t')
    print ""
else:
    for alt in alternatives:
        sys.stdout.write(getAlt(alt) + '_AVG\t' + getAlt(alt) + '_STDDEV\t')
    print ""

for bench in benchmarks:
    sys.stdout.write(bench + '\t')
    for contract in contracts:
        if not 'ReconfigurationTime' in args.field:
            for alt in alternatives:
                printField(bench, contract, alt, args.field)
            print ""
        else:
            for alt in alternativesReconfTime:
                printField(bench, contract, alt, args.field)
            print ""

