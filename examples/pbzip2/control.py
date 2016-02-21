import time
from subprocess import check_output
import math
import subprocess
import shlex
import socket
import os
import os.path
import subprocess, signal

def kill():
    p = subprocess.Popen(['ps', '-A'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.splitlines():
        if 'pbzip2' in line:
            pid = int(line.split(None, 1)[0])
            os.kill(pid, signal.SIGKILL)
            print "Killing"

def hang():
    if os.path.isfile("summary.csv"):
        run = "wc -l summary.csv"
        process = subprocess.Popen(shlex.split(run), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = process.communicate()
        return int(out.split(' ')[0]) == 2
    else:
        return 0

def get_pid(name):
    return check_output(["pidof",name]).split("\n")[0]

while 1:
    time.sleep(60)
    if hang():
        time.sleep(10)
        if hang():
            kill()
