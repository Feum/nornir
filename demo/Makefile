CXX                  = g++
CXXFLAGS             = -Wall -g -pedantic -fpermissive
OPTIMIZE_FLAGS       = -O3 -finline-functions 
INCS                 = -I../../Mammut
LDFLAGS              = -L../../Mammut/mammut
LDLIBS               = -lm -pthread -lrt -lmammut -lprotobuf-lite
TARGET               = demo-farm demo-parfor

.PHONY: clean cleanall

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) --std=c++11
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

