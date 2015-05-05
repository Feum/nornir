#!/bin/bash
NUMCORES=22
FREQUENCIES=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies | tr ' ' '\n' | sort | tr '\n' ' ')

for (( i=0; i<=47; i++ ))
do
    cpufreq-set -g userspace -c $i
done

rm -rf results.txt
echo "Workers Frequency JoulesCpu JoulesCores CompletionTime" | tee results.txt

for (( i=2; i<=$NUMCORES; i++ ))
do
  for j in $FREQUENCIES
  do
    for (( k=0; k<=47; k++ ))
     do
      cpufreq-set -c $k -f $j
    done
    CT=$(./pbzip2_ff /home/desensi/enwiki-latest-pages-articles.xml -f -k -p$i | grep 'CompletionTimeSecs' | cut -d ':' -f 2)
    JOULESCPU=$(tail -1 energy.txt | cut -d ' ' -f 1)
    JOULESCORES=$(tail -1 energy.txt | cut -d ' ' -f 2)
    echo $i $j $JOULESCPU $JOULESCORES $CT | tee --append results.txt
    rm sample.pcap.bz2
    mv stats.txt stats$i$j.txt
    mv energy.txt energy$i$j.txt
  done
done

for (( i=0; i<=47; i++ ))
do
    cpufreq-set -g performance -c $i
done