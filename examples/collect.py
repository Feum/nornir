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
fields['CalibrationWattsAvg'] = 22
fields['CalibrationWattsStddev'] = 23
fields['ReconfigurationTimeWorkersAvg'] = 24
fields['ReconfigurationTimeWorkersStddev'] = 25
fields['ReconfigurationTimeFrequencyAvg'] = 26
fields['ReconfigurationTimeFrequencyStddev'] = 27
fields['ReconfigurationTimeTotalAvg'] = 28
fields['ReconfigurationTimeTotalStddev'] = 29

alternatives = ('LIMARTINEZ_SLEEP_SMALL', 'REGRESSION_LINEAR_RANDOM_FAST_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_SLEEP_LATENCY', 'REGRESSION_LINEAR_HALTON_FAST_SPINNING', 'REGRESSION_LINEAR_HALTON_FAST_MAPPING_SLEEP_SMALL', 'INTEL_SLEEP_SMALL')
alternativesReconfTime = ('REGRESSION_LINEAR_HALTON', 'REGRESSION_LINEAR_HALTON_FAST')
alternativesMandel = ('LIMARTINEZ', 'REGRESSION_LINEAR_HALTON_FAST_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE10_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE20_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_CONSERVATIVE30_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_AGING4_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_AGING8_SLEEP_SMALL', 'REGRESSION_LINEAR_HALTON_FAST_AGING12_SLEEP_SMALL')

def printField(bench, contract, alt, field):
    try:
        fh = open(bench + "/" + contract + "_" + alt + "/results.csv")
        fh.readline() # Skip header
        values = []
        variances = []

        for i in xrange(10, 100, 20):
            fieldsLine = fh.readline()
            fieldValue = float(fieldsLine.split('\t')[fields[field]])
            if 'CalibrationWatts' in field and 'POWER_BUDGET' in contract:
                budget = float(fieldsLine.split('\t')[fields['PrimaryRequired']])
                fieldValue = ((fieldValue - budget)/budget)*100.0
                
            if (field != 'PrimaryLossAvg' and field != 'SecondaryLossAvg') or fieldValue != 0:
                values.append(fieldValue)
                variances.append(float(fieldsLine.split('\t')[fields[field] + 1]) ** 2) # ^2 because in the field there is stddev e we need the variance

        if field == 'PrimaryLossCnt':
            sys.stdout.write(str(np.sum(values)) + '\t0\t')
        else:
            if len(values):
                avg = np.average(values)
                std = math.sqrt(np.sum(variances)/len(variances)) # Pooled variance formula
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
benchmarks = ('blackscholes', 'canneal', 'pbzip2', 'videoprocessing', 'simple_mandelbrot')
contracts = ('PERF_BANDWIDTH', 'POWER_BUDGET')
fields 

parser = argparse.ArgumentParser(description='Runs the application to check the accuracy of the reconfiguration.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-b', '--benchmark', help='Benchmark.', required=False)
parser.add_argument('-c', '--contract', help='Contract type.', required=True)
parser.add_argument('-f', '--field', help='Field type (PrimaryLossCnt, PrimaryLossAvg, SecondaryLossAvg, CalibrationTimePercAvg, ReconfigurationTimeWorkersAvg).', required=True)
args = parser.parse_args()


if args.benchmark is not None:
    benchmarks = ([args.benchmark])
    if 'simple_mandelbrot' in args.benchmark:
        alternatives = alternativesMandel

if args.contract is not None:
    contracts = ([args.contract])

sys.stdout.write('#Bench\t')

for alt in alternatives:
    sys.stdout.write(getAlt(alt) + '_AVG\t' + getAlt(alt) + '_STDDEV\t')
print ""

for bench in benchmarks:
    if 'mandelbrot' in bench:
        sys.stdout.write('mandelbrot\t')
    elif 'videoprocessing' in bench:
        sys.stdout.write('denoiser\t')
    else:
        sys.stdout.write(bench + '\t')
    for contract in contracts:
        for alt in alternatives:
            printField(bench, contract, alt, args.field)
        print ""


