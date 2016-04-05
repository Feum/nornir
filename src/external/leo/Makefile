.PHONY = all

all: leo.o

leo.o: leo.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 -c leo.cpp
clean: 
	-rm -fr *.o *~
cleanall:
	make clean
