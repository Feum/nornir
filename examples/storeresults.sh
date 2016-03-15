#!/bin/bash

mkdir -p RESULTS
for bench in videoprocessing blackscholes canneal simple_mandelbrot pbzip2
do
    mkdir -p RESULTS/$bench
    cp -R $bench/PERF* RESULTS/$bench/ 
    cp -R $bench/POWER* RESULTS/$bench/
done
