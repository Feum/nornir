import math
import subprocess
import shlex
import socket
import os
import shutil 
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

#fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
fractions = [0.9, 1.0]

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
    parametersFile.write("<persistenceValue>3</persistenceValue>\n")
    parametersFile.write("<strategySmoothing>MOVING_AVERAGE</strategySmoothing>\n")
    parametersFile.write("<contractType>" + contractType + "</contractType>\n")
    parametersFile.write("<" + fieldName + ">" + str(fieldValue) + "</" + fieldName + ">\n")
    parametersFile.write("<strategyPolling>SLEEP_SMALL</strategyPolling>\n")
    parametersFile.write("</adaptivityParameters>\n")
    parametersFile.close()

    run = "./pbzip2_ff -f -k -p24 /home/desensi/enwiki-20151201-abstract.xml"
    process = subprocess.Popen(shlex.split(run), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = process.communicate()

    print out 
    print err

    result = getLastLine("summary.csv")

    dir_results = contractType
    for outfile in ["calibration.csv", "stats.csv", "summary.csv", "parameters.xml"]:
        os.rename(outfile, contractType + "/" + str(fraction) + "." + outfile)

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
    fh = open("REPARA_results.csv", "r")
    bestTime = 9999999999
    for line in fh:
        #Workers Frequency Time Watts
        if line[0] != '#':
            print line
            fields = line.split("\t")
            curtime = float(fields[2])
            curpower = float(fields[3])
            if curtime < bestTime and curpower <= float(power):
                bestTime = curtime
    fh.close()
    return bestTime

############################################################################################################

shutil.rmtree('PERF_COMPLETION_TIME', ignore_errors=True)
shutil.rmtree('POWER_BUDGET', ignore_errors=True)
os.makedirs('PERF_COMPLETION_TIME')
os.makedirs('POWER_BUDGET')

perfFile = open("PERF_COMPLETION_TIME/results.csv", "w")
powerFile = open("POWER_BUDGET/results.csv", "w")
perfFile.write("#Fraction\tPrimaryRequired\tPrimaryFound\tSecondaryOptimal\tSecondaryFound\tLoss\n")
powerFile.write("#Fraction\tPrimaryRequired\tPrimaryFound\tSecondaryOptimal\tSecondaryFound\tLoss\n")

for f in fractions:
    targetTime = timeMin + timeRange*f
    ct, watts = run("PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, f)
    opt = getOptimalTimeBound(targetTime)
    loss = ((watts - opt) / opt) * 100.0
    perfFile.write(str(f) + "\t" + str(targetTime) + "\t" + str(ct)  + "\t" + str(opt) + "\t" + str(watts) + "\t" + str(loss) + "\n")

    targetPower = powerMin + powerRange*f
    ct, watts = run("POWER_BUDGET", "powerBudget", targetPower, f)
    opt = getOptimalPowerBound(targetPower)
    loss = ((ct - opt) / opt) * 100.0
    powerFile.write(str(f) + "\t" + str(targetPower) + "\t" + str(watts) + "\t" + str(opt) + "\t" + str(ct) + "\t" + str(loss) + "\n")

perfFile.close()
powerFile.close()
