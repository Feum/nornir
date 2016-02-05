#!/bin/bash
while kill -0 24384; do
    sleep 60
done


for bench in simple_mandelbrot canneal #pbzip2 blackscholes #videoprocessing #simple_mandelbrot canneal
do
    (cd $bench; python control.py) &
    CPID=$!
    cd $bench; python runreconf.py --prediction REGRESSION_LINEAR --contract PERF_COMPLETION_TIME; cd ..
    cd $bench; python runreconf.py --prediction REGRESSION_LINEAR --contract POWER_BUDGET; cd ..
    cd $bench; python runreconf.py --prediction LIMARTINEZ --contract PERF_COMPLETION_TIME; cd ..
    kill -9 $CPID
done
