#!/bin/bash
echo "Counting the number of frames of video: " $1
#FRAMES=$(ffmpeg -i "$1" -vcodec copy -acodec copy -f null /dev/null 2>&1 | grep 'frame=' | cut -f 2 -d '=' | cut -f 1 -d ' ')
FRAMES=$(ffprobe -v error -count_frames -select_streams v:0 -show_entries stream=nb_read_frames -of default=nokey=1:noprint_wrappers=1 "$1")
echo $FRAMES
