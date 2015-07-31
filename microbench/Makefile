CXX                  = g++
CXXFLAGS             = -Wall -g -pedantic -fpermissive 
OPTIMIZE_FLAGS       = -O3 -finline-functions 
INCS                 = 
LDFLAGS              = 
LDLIBS               = -pthread -lrt
TARGET               = ticksPerNs

.PHONY: clean cleanall

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) --std=c++11
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

