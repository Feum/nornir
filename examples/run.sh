#!/bin/bash

#while kill -0 47319; do
#    sleep 60
#done

declare -a BENCHS=("videoprocessing" "canneal" "blackscholes" "simple_mandelbrot" "pbzip2")

for bench in "${BENCHS[@]}"
do
    (cd $bench; python control.py) &
    CPID=$!
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_BANDWIDTH -f
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f
done

for bench in "${BENCHS[@]}"
do  
    #(cd $bench; python control.py) &
    #CPID=$!
    ##python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -k RANDOM
    #python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SLEEP_LATENCY
    #python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -l SPINNING

    python runreconf.py -b $bench -p LIMARTINEZ -c PERF_BANDWIDTH
done

for bench in "${BENCHS[@]}"
do
    ##python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -k RANDOM
    #python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SLEEP_LATENCY
    #python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_COMPLETION_TIME -f -l SPINNING

    python runreconf.py -b $bench -p REGRESSION_LINEAR -c POWER_BUDGET -f -w MAPPING
    python runreconf.py -b $bench -p REGRESSION_LINEAR -c PERF_BANDWIDTH -f -w MAPPING
done

for bench in "${BENCHS[@]}"
do
    python runreconf.py -b $bench -p INTEL -c POWER_BUDGET
done


