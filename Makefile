export NORNIR_PATH           = /usr/local
export NORNIR_PATH_LIB       = $(NORNIR_PATH)/lib
export NORNIR_PATH_INCLUDE   = $(NORNIR_PATH)/include/nornir

export CC                    = gcc
export CXX                   = g++ 
export OPTIMIZE_FLAGS        = -finline-functions  -O3 #-DPOOL
export DEBUG_FLAGS           = #-DDEBUG_PREDICTORS -DDEBUG_SELECTORS #-DDEBUG_DF_STREAM -DDEBUG_NODE -DDEBUG_KNOB -DDEBUG_MANAGER
export CXXFLAGS              = -Wall -pedantic --std=c++11 $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)
export LDLIBS                =  -lnornir -pthread -lrt -lm -lmlpack -llapack -lblas -lgsl -lgslcblas -larmadillo -lorlog -lanl
export INCS                  = -I$(realpath ./src/external/fastflow) -I/usr/include/libxml2
export LDFLAGS               = -L$(realpath .)/src -L$(realpath .)/src/external/orlog/src

.PHONY: all demo clean cleanall install uninstall microbench bin

all:
	python submodules_init.py
	git submodule foreach git pull -q origin master
	$(MAKE) -C src
	$(MAKE) -C microbench
	$(MAKE) -C microbench checksupported
clean: 
	$(MAKE) -C src clean
	$(MAKE) -C demo clean
	$(MAKE) -C examples clean
	$(MAKE) -C microbench clean
	$(MAKE) -C bin clean
demo:
	$(MAKE) -C demo 
	$(MAKE) -C examples
bin:
	$(MAKE) -C bin
cleanall:
	$(MAKE) -C src cleanall
	$(MAKE) -C demo cleanall
	$(MAKE) -C examples cleanall
	$(MAKE) -C microbench cleanall
	$(MAKE) -C bin cleanall
install:
	$(MAKE) -C src install
uninstall:
	$(MAKE) -C src uninstall
microbench:
	$(MAKE) -C microbench microbench
