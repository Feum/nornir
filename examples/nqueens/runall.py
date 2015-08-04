#!/usr/bin/env python 
from subprocess import call
import multiprocessing, os, errno, shlex

output='results.csv'

def silentRemove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def setFrequency(f):
    for coreId in xrange(0, multiprocessing.cpu_count()):
        call(shlex.split("cpufreq-set -g userspace -c " + str(coreId))) 
        call(shlex.split("cpufreq-set -f " + f + " -c " + str(coreId)))
    return

with open("/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies") as f:
    frequencies = f.read().split()
frequencies.reverse()
# Remove turboboost frequency
if frequencies[-1][3] == '1':
    del frequencies[-1]

silentRemove(output)
with open(output, 'a') as results:
    results.write('Workers\tFrequency\tWattsCores\tBandwidth\tCompletionTime\n')

for numWorkers in xrange(3,23):
    for frequency in frequencies:
        setFrequency(frequency)
        call(["./nq_ff", "18 " + str(numWorkers) + " 10"])
        with open("summary.csv") as summary:
            fields = summary.readlines()[-1].split('\t')
        with open(output, 'a') as results:
            results.write(str(numWorkers)+"\t"+frequency+"\t"+
                          fields[0]+"\t"+fields[1]+"\t"+fields[2])
            results.write("\n");
            
        

