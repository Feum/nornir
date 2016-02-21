#!/bin/bash

for bench in canneal pbzip2
do
    (cd $bench; python control.py) &
    CPID=$!
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -k RANDOM
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -k HALTON
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SLEEP_LATENCY
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SPINNING

    python runreconf.py -b $bench -p LIMARTINEZ -c PERF_COMPLETION_TIME

    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -k RANDOM
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -k HALTON
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SLEEP_LATENCY
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SPINNING
done 

for bench in videoprocessing blackscholes canneal pbzip2 simple_mandelbrot
do  
    (cd $bench; python control.py) &
    CPID=$!
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -w MAPPING
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -w MAPPING
done

: '
for bench in videoprocessing blackscholes canneal pbzip2 simple_mandelbrot
do  
    (cd $bench; python control.py) &
    CPID=$!
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -k RANDOM
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -k HALTON
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SLEEP_LATENCY
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SPINNING

    python runreconf.py -b $bench -p LIMARTINEZ -c PERF_COMPLETION_TIME

    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -k RANDOM
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -k HALTON
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SLEEP_LATENCY
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SPINNING
done
'

: '
for bench in videoprocessing #simple_mandelbrot #videoprocessing blackscholes canneal pbzip2
do
    (cd $bench; python control.py) &
    CPID=$!

    python runreconf.py -b $bench -p LIMARTINEZ -c PERF_BANDWIDTH 
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f

    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f 

#    python runreconf.py -b $bench -p LIMARTINEZ -c PERF_BANDWIDTH -r REPARA_results_phases.csv
#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_BANDWIDTH -r REPARA_results_phases.csv
#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_BANDWIDTH -f -r REPARA_results_phases.csv

#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -r REPARA_results_phases.csv
#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -r REPARA_results_phases.csv
    kill -9 $CPID
done

for bench in blackscholes canneal pbzip2 #videoprocessing simple_mandelbrot
do
    (cd $bench; python control.py) &
    CPID=$!
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SLEEP_LATENCY
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SLEEP_LATENCY

    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SPINNING
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SPINNING
    kill -9 $CPID
done

#for bench in simple_mandelbrot #videoprocessing blackscholes canneal pbzip2
#do
#    (cd $bench; python control.py) &
#    CPID=$!
#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -o 10
#    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -o 30
#    kill -9 $CPID
#done

'
