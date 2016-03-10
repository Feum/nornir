#!/bin/sh
export LD_PRELOAD=/home/torquati/Hoard/src/libhoard.so
./power_gov -r POWER_LIMIT -s 2000 -d PP0 
./power_gov -r CLAMPING_LIMIT -s 0 -d PP0 
./power_gov -r ENABLE_POWER_LIMIT -s 0 -d PP0

