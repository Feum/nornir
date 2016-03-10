#!/bin/sh

export LD_PRELOAD=/home/torquati/Hoard/src/libhoard.so

PER_SOCKET_CAP=$(echo "$1" | awk '{printf "%.2f \n", $1/2}')
echo "Cap per socket: " $PER_SOCKET_CAP

./power_gov -r ENABLE_POWER_LIMIT -s 1 -d PP0 
./power_gov -r POWER_LIMIT -s $PER_SOCKET_CAP -d PP0 
./power_gov -r TIME_WINDOW -s 1 -d PP0 
./power_gov -r CLAMPING_LIMIT -s 1 -d PP0 
