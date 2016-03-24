#!/bin/bash

rm -rf archdata.xml
# Header
echo "<?xml version="1.0" encoding="UTF-8"?>" >> archdata.xml
echo "<archData>" >> archdata.xml

# TickPerNs
echo -n "    " >> archdata.xml
./ticksPerNs >> archdata.xml

# Voltage table
if [ ! -f voltageTable.txt ]; then
	if [ "$(id -u)" != "0" ]; then
		echo "[[[[ ATTENTION ]]]] You need to run this command with sudo (this is needed in order to compute the voltage table)."
		exit 1
	fi
	./voltageTable
fi
echo -n "    " >> archdata.xml
echo "<voltageTableFile>voltage.txt</voltageTableFile>" >> archdata.xml

# Footer
echo "</archData>" >> archdata.xml

# Copy all
for i in ../examples/*; do 
	cp voltage.txt $i; 
	cp archdata.xml $i;
done;

cp voltage.txt ../demo
cp archdata.xml ../demo
