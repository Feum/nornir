LDFLAGS              = -L/home/daniele/Mammut/mammut
LDLIBS               = -lm -pthread -lrt -lmammut -lprotobuf-lite
TARGET               = demo-fastflow

.PHONY: clean cleanall

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

