LDLIBS                 += -lpcap
TARGET               = ffProbe

.PHONY: all clean cleanall install uninstall
.SUFFIXES: .cpp .o

all: $(TARGET)

ffProbe.o: ffProbe.cpp *.hpp
	$(CXX) $(INCS) $(CXXFLAGS) $(OPTIMIZE_FLAGS) -c -o $@ $<

ffProbe: ffProbe.o
	$(CXX) $< -o $@ $(LDFLAGS) $(LDLIBS)
	
clean: 
	-rm -fr *.o *~
	
cleanall: clean
	-rm -fr $(TARGET)
	
install:
	cp ./ffProbe /usr/local/bin/ffProbe
	
uninstall:
	rm -fr /usr/local/bin/ffProbe
