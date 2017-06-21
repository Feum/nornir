TARGET = monitor.o manager-multi.o manager.o knob.o trigger.o node.o parameters.o selectors.o predictors.o explorers.o configuration.o stats.o

.PHONY: all clean cleanall install uninstall lib

all:
	$(MAKE) -C external
	$(MAKE) -C dataflow
	$(MAKE) lib
%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(OPTIMIZE_FLAGS) -c -o $@ $<
clean: 
	-rm -fr *.o *~
	$(MAKE) -C dataflow clean
	$(MAKE) -C external clean
cleanall:
	-rm -fr *.o *~ *.a
	$(MAKE) -C dataflow cleanall
#	$(MAKE) -C external cleanall
lib: $(TARGET)
	rm -rf libnornir.a
	ar rc libnornir.a $(TARGET) dataflow/*.o
	ranlib libnornir.a
	ar -M < buildlib.ar
install:
	mkdir -p $(NORNIR_PATH_INCLUDE)
	find . -name '*.h' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	find . -name '*.hpp' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	find . -name '*.tpp' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	cp ./libnornir.a $(NORNIR_PATH_LIB)/libnornir.a
uninstall:
	rm -r $(NORNIR_PATH_INCLUDE)/*
	rm $(NORNIR_PATH_LIB)/libnornir.a

