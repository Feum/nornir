INCS                 += -I$(realpath ../src/external/fastflow)
TARGET               = ticksPerNs

.PHONY: all clean cleanall microbench

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) 
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)
microbench:
	cp ../src/external/Mammut/demo/energy/voltageTable.cpp .
	$(CXX) $(CXXFLAGS) -o voltageTable voltageTable.cpp $(INCS) -I../src/external/Mammut/  -L../src/external/Mammut/mammut $(LDFLAGS) $(LDLIBS)
	@echo ""
	@echo ""
	@echo "=============================== ATTENTION ==============================="
	@echo "| Sudo rights are needed to compute the voltage table. If voltage table |" 
	@echo "| is not computed, it will not be possible to predict power consumption |"
	@echo "| of the applications.                                                  |" 
	@echo "========================================================================="
	sudo sh runmicrobenchs.sh
	@echo "=============================== ATTENTION ==============================="
	@echo "| archdata.xml file has been produced in the microbench folder. Copy    |" 
	@echo "| this file to the folders of the applications you want to run.         |"
	@echo "========================================================================="
