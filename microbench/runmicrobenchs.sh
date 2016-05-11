#!/bin/bash

# When a new parameter is added to the archdata.xml configuration file.
# 1. Change the CURRENT_VERSION_VARIABLE
# 2. Insert another condition in the if for checking if the new field
#    has been correctly created (e.g. 'grep "ticksPerNs" ... ')

CURRENT_VERSION="1.0.0"

CONFPATH_ROOT=$XDG_CONFIG_HOME
if [ -z "$VAR" ];
then
	CONFPATH_ROOT=$HOME"/.config"
fi

echo "Creating "$CONFPATH_ROOT"/nornir folder."
mkdir -p $CONFPATH_ROOT"/nornir"

CONFPATH_FILE=$CONFPATH_ROOT"/nornir/archdata.xml"
CONFPATH_VOLTAGE=$CONFPATH_ROOT"/nornir/voltage.csv"
CONFPATH_VERFILE=$CONFPATH_ROOT"/nornir/version.csv"

if grep $CURRENT_VERSION "$CONFPATH_VERFILE" >/dev/null 2>&1 && \
   grep "ticksPerNs" "$CONFPATH_FILE" >/dev/null 2>&1 && \
   [ -f $CONFPATH_VOLTAGE ] && \
   [ -f $CONFPATH_FILE ]; 
then
	echo "Configuration file for version "$CURRENT_VERSION " already created."
else
	echo $CURRENT_VERSION > $CONFPATH_VERFILE
	rm -rf $CONFPATH_FILE
	# Header
	echo "<?xml version="1.0" encoding="UTF-8"?>" >> $CONFPATH_FILE
	echo "<archData>" >> $CONFPATH_FILE
	
	# TickPerNs
	echo -n "    " >> $CONFPATH_FILE
	./ticksPerNs >> $CONFPATH_FILE
	
	# Voltage table
	if [ ! -f $CONFPATH_VOLTAGE ]; then
		if [ ! -f voltageTable.txt ]; then
			if [ "$(id -u)" != "0" ]; then
				echo "[[[[ ATTENTION ]]]] You need to run this command with sudo (this is needed in order to compute the voltage table)."
				exit 1
			fi
			./voltageTable
		fi
		cp voltageTable.txt $CONFPATH_VOLTAGE
	fi
	
	# Footer
	echo "</archData>" >> $CONFPATH_FILE
fi