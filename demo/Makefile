TARGET               = farm farm_ffcompat

.PHONY: clean cleanall

all: $(TARGET)

INCS += -I$(..)

%: %.cpp ../src/libnornir.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

