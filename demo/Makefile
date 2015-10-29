CXX                  = g++
CXXFLAGS             = -Wall -g -pedantic -fpermissive -DDEBUG_FARM -DDEBUG_NODE
OPTIMIZE_FLAGS       = -O3 -finline-functions 
INCS                 = -I../../Mammut -I/usr/include/libxml2
LDFLAGS              = -L../../Mammut/mammut
LDLIBS               = -lm -pthread -lrt -lmammut -lprotobuf-lite -lmlpack -llapack -lblas -lgsl -lgslcblas
TARGET               = demo-farm demo-farm-sendout

.PHONY: clean cleanall

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) --std=c++11
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

