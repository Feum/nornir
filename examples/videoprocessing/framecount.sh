#!/bin/bash
echo "Counting the number of frames of video: " $1
FRAMES=$(ffmpeg -i "$1" -vcodec copy -acodec copy -f null /dev/null 2>&1 | grep 'frame=' | cut -f 2 -d ' ')
echo $FRAMES
