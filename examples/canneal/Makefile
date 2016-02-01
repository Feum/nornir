# Makefile for parallel simulated annealer

TARGET = canneal
LDLIBS += $(LIBS) -lm

CXXFLAGS += -DENABLE_FF -DENABLE_THREADS -O3

all: 
	$(CXX) $(CXXFLAGS) annealer_thread.cpp -c -o annealer_thread.o
	$(CXX) $(CXXFLAGS) rng.cpp -c -o rng.o
	$(CXX) $(CXXFLAGS) netlist.cpp -c -o netlist.o
	$(CXX) $(CXXFLAGS) $(INCS) main.cpp -c -o main.o
	$(CXX) $(CXXFLAGS) netlist_elem.cpp -c -o netlist_elem.o
	$(CXX) $(CXXFLAGS)  $(LDFLAGS) *.o $(LDLIBS) -o $(TARGET)

clean:
	rm -f *.o $(TARGET)


