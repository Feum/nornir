import math
import subprocess
import shlex
import socket
import os
import shutil
import argparse 
from subprocess import Popen, PIPE, STDOUT

RESULTS_FILE = 'REPARA_results.csv'

cmd = 'cat ' + RESULTS_FILE + ' | tail -n +2 | cut -f 3 | sort -n | head -1'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
timeMin = float(p.stdout.read())

cmd = 'cat ' + RESULTS_FILE + ' | tail -n +2 | cut -f 3 | sort -n | tail -1'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
timeMax = float(p.stdout.read())

cmd = 'cat ' + RESULTS_FILE + ' | tail -n +2 | cut -f 4 | sort -n | head -1'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
powerMin = float(p.stdout.read())

cmd = 'cat ' + RESULTS_FILE + ' | tail -n +2 | cut -f 4 | sort -n | tail -1'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
powerMax = float(p.stdout.read())

timeRange = timeMax - timeMin
powerRange = powerMax - powerMin

print "TimeMin: " + str(timeMin) + " TimeMax: " + str(timeMax)
print "PowerMin: " + str(powerMin) + " PowerMax: " + str(powerMax)

fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

def getLastLine(fileName):
    fh = open(fileName, "r")
    for line in fh:
        pass
    last = line
    fh.close()
    return last


def run(contractType, fieldName, fieldValue, fraction):
    parametersFile = open("parameters.xml", "w")
    parametersFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    parametersFile.write("<adaptivityParameters>\n")
    parametersFile.write("<qSize>4</qSize>\n")
    parametersFile.write("<samplingInterval>500</samplingInterval>\n")

    parametersFile.write("<strategyPrediction>" + args.prediction + "</strategyPrediction>\n")

    parametersFile.write("<strategyPersistence>TASKS</strategyPersistence>\n")
    parametersFile.write("<persistenceValue>1000</persistenceValue>\n")
    parametersFile.write("<strategySmoothing>MOVING_AVERAGE</strategySmoothing>\n")

    parametersFile.write("<contractType>" + contractType + "</contractType>\n")
    parametersFile.write("<" + fieldName + ">" + str(fieldValue) + "</" + fieldName + ">\n")
    parametersFile.write("<strategyPolling>SLEEP_SMALL</strategyPolling>\n")
    parametersFile.write("</adaptivityParameters>\n")
    parametersFile.close()

    run = "./blackscholes 23 inputs/in_10M.txt tmp.txt"
    process = subprocess.Popen(shlex.split(run), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = process.communicate()

    result = getLastLine("summary.csv")

    dir_results = "nodir"
    if contractType == 'PERF_COMPLETION_TIME':
        dir_results = perfdir
    elif contractType == 'POWER_BUDGET':
        dir_results = powerdir

    for outfile in ["calibration.csv", "stats.csv", "summary.csv", "parameters.xml"]:
        os.rename(outfile, dir_results + "/" + str(fraction) + "." + outfile)

    # Returns completion time and watts
    return float(result.split('\t')[2]), float(result.split('\t')[0])

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

if args.contract == "PERF_COMPLETION_TIME":
    perfdir = 'PERF_COMPLETION_TIME_' + args.prediction
    shutil.rmtree(perfdir, ignore_errors=True)
    os.makedirs(perfdir)
    perfFile = open(perfdir + "/results.csv", "w")
    perfFile.write("#Fraction\tPrimaryRequired\tPrimaryFound\tSecondaryOptimal\tSecondaryFound\tLoss\n")
elif args.contract == "POWER_BUDGET":
    powerdir = 'POWER_BUDGET_' + args.prediction
    shutil.rmtree(powerdir, ignore_errors=True)
    os.makedirs(powerdir)
    powerFile = open(powerdir + "/results.csv", "w")
    powerFile.write("#Fraction\tPrimaryRequired\tPrimaryFound\tSecondaryOptimal\tSecondaryFound\tLoss\n")

for f in fractions:
    if args.contract == "PERF_COMPLETION_TIME":
        targetTime = timeMin + timeRange*f
        ct, watts = run("PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, f)
        opt = getOptimalTimeBound(targetTime)
        loss = ((watts - opt) / opt) * 100.0
        perfFile.write(str(f) + "\t" + str(targetTime) + "\t" + str(ct)  + "\t" + str(opt) + "\t" + str(watts) + "\t" + str(loss) + "\n")
        perfFile.flush()
        os.fsync(perfFile.fileno())
    elif args.contract == "POWER_BUDGET":
        targetPower = powerMin + powerRange*f
        ct, watts = run("POWER_BUDGET", "powerBudget", targetPower, f)
        opt = getOptimalPowerBound(targetPower)
        loss = ((ct - opt) / opt) * 100.0
        powerFile.write(str(f) + "\t" + str(targetPower) + "\t" + str(watts) + "\t" + str(opt) + "\t" + str(ct) + "\t" + str(loss) + "\n")
        powerFile.flush()
        os.fsync(powerFile.fileno())


if args.contract == "PERF_COMPLETION_TIME":
    perfFile.close()
elif args.contract == "POWER_BUDGET":
    powerFile.close()
