MAMMUT_ROOT                 = $(realpath ./src/external/Mammut)

export GRAPE_PATH           = /usr/local
export GRAPE_PATH_LIB       = $(GRAPE_PATH)/lib
export GRAPE_PATH_INCLUDE   = $(GRAPE_PATH)/include/grape

export CC                    = gcc
export CXX                   = g++ 
export OPTIMIZE_FLAGS        = -finline-functions -O3 #-O0 -g 
export DEBUG_FLAGS           = #-DDEBUG_NODE -DDEBUG_KNOB -DDEBUG_PREDICTORS -DDEBUG_MANAGER -DDEBUG_NODE
export CXXFLAGS              = -Wall -pedantic --std=c++11 -DFF_TASK_CALLBACK -DTRACE_FASTFLOW $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)
export LDLIBS                =  -lgrape -pthread -lrt -lm -lmlpack -llapack -lblas -lgsl -lgslcblas 
export INCS                  = -I$(realpath ./src/external/fastflow) -I$(MAMMUT_ROOT) -I/usr/include/libxml2
export LDFLAGS               = -L$(MAMMUT_ROOT)/mammut -L$(realpath .)/src

.PHONY: all demo clean cleanall install uninstall

all:
	python submodules_init.py
	git submodule foreach git pull -q origin master
	$(MAKE) -C src
clean: 
	$(MAKE) -C src clean
demo:
	$(MAKE) -C demo 
	$(MAKE) -C examples
cleanall:
	$(MAKE) -C src cleanall
	$(MAKE) -C demo cleanall
	$(MAKE) -C examples cleanall
install:
	$(MAKE) -C src install
uninstall:
	$(MAKE) -C src uninstall
