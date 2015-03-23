export INCS                  = -I/home/daniele/Mammut
export CC                    = gcc
export CXX                   = g++
export OPTIMIZE_FLAGS        = -O3 -finline-functions 
export CXXFLAGS              = -Wall -g -pedantic 

.PHONY: all

all: 
	make -C demo 
