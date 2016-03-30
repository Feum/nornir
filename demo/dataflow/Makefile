TARGET               = asciiMandelbrot

.PHONY: clean cleanall

all: $(TARGET)

INCS += -I$(realpath ../../src)

%: %.cpp ../../src/libnornir.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

