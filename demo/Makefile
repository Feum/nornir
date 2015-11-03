TARGET               = demo-farm demo-farm-sendout

.PHONY: clean cleanall

all: $(TARGET)

%: %.cpp ../src/libgrape.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

