#!/bin/bash
# Converts all the .svg files in the folder to .png, changing the color.
# $1: Destination color (without #)

echo "Creating with color " $1
for i in *.svg
do	
	sed -e "s/#000000/#$1/" $i > $i.color
	inkscape --export-png="${i%.*}".$1.png --export-dpi=1200 --export-background-opacity=0 --without-gui $i.color
	rm $i.color
done
