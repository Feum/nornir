import math
import subprocess
import shlex
import socket
import os
import shutil 

timeMin = 72.677
timeMax = 1564.75
timeRange = timeMax - timeMin

powerMin = 24.9002
powerMax = 110.6
powerRange = powerMax - powerMin

#fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

fractions = [0.1]

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
    parametersFile.write("<samplingInterval>1000</samplingInterval>\n")
    parametersFile.write("<persistenceValue>6</persistenceValue>\n")
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
    fh = open("REPARA_results.csv", "r")
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
                bestPower = curpower
    fh.close()
    return bestTime

############################################################################################################

shutil.rmtree('PERF_COMPLETION_TIME', ignore_errors=True)
shutil.rmtree('POWER_BUDGET', ignore_errors=True)
os.makedirs('PERF_COMPLETION_TIME')
os.makedirs('POWER_BUDGET')

perfFile = open("PERF_COMPLETION_TIME/results.csv", "w")
powerFile = open("POWER_BUDGET/results.csv", "w")
perfFile.write("#Fraction\tSolution\tLoss")
powerFile.write("#Fraction\tSolution\tLoss")

for f in fractions:
    targetTime = timeMin + timeRange*f
    ct, watts = run("PERF_COMPLETION_TIME", "requiredCompletionTime", targetTime, f)
    opt = getOptimalTimeBound(ct)
    loss = ((watts - opt) / opt) * 100.0
    perfFile.write(str(f) + str(watts) + str(loss) + "\n")

    targetPower = powerMin + powerRange*f
    ct, watts = run("POWER_BUDGET", "powerBudget", targetPower, f)
    opt = getOptimalPowerBound(watts)
    loss = ((ct - opt) / opt) * 100.0
    powerFile.write(str(f) + str(ct) + str(loss) + "\n")

perfFile.close()
powerFile.close()
