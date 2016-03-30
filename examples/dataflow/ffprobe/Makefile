LDLIBS               += -lpcap
TARGET               = ffProbe
CXXFLAGS             += -DPOOL

.PHONY: all clean cleanall install uninstall
.SUFFIXES: .cpp .o

all: $(TARGET)

ffProbe.o: ffProbe.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<

ffProbe: ffProbe.o
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean: 
	-rm -fr *.o *~

cleanall: clean
	-rm -fr $(TARGET)

install:
	cp ./ffProbe /usr/local/bin/ffProbe

uninstall:
	rm -fr /usr/local/bin/ffProbe
