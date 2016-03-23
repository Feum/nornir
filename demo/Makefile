TARGET               = farm farm_ffcompat farm_accelerator_1 farm_accelerator_2 farm_accelerator_3 farm_accelerator_4

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

