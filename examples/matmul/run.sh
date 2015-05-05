#!/bin/bash
NUMCORES=22
FREQUENCIES=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies | tr ' ' '\n' | sort | tr '\n' ' ')

for (( i=0; i<=47; i++ ))
do
    cpufreq-set -g userspace -c $i
done

rm -rf results.txt
echo "Workers Frequency WattsCpu WattsCores WattsOffCores CompletionTime" | tee results.txt

for (( i=2; i<=$NUMCORES; i++ ))
do
for j in $FREQUENCIES
do
    for (( k=0; k<=47; k++ ))
     do
      cpufreq-set -c $k -f $j
    done
    CT=$(./ff_streamMM 512 2000 $i | grep "Total time" | cut -d ' ' -f 3 | cut -d ' ' -f 1)
    WATTSCPU=$(tail -1 energy.txt | cut -d ' ' -f 1)
    WATTSCORES=$(tail -1 energy.txt | cut -d ' ' -f 2)

    WATTSOFFCORES=$(awk "BEGIN {printf \"%.2f\",${WATTSCPU}-${WATTSCORES}}")
    echo $i $j $WATTSCPU $WATTSCORES $WATTSOFFCORES $CT | tee --append results.txt
    mv stats.txt stats$i$j.txt
    mv energy.txt energy$i$j.txt
  done
done

for (( i=0; i<=47; i++ ))
do
    cpufreq-set -g performance -c $i
done