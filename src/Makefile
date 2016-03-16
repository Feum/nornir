TARGET = knob.o trigger.o node.o parameters.o predictors.o configuration.o utils.o

.PHONY: all clean cleanall install uninstall lib

all:
	$(MAKE) -C external
	$(MAKE) lib
%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(OPTIMIZE_FLAGS) -c -o $@ $<
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ *.a
lib: $(TARGET)
	ar rc libnornir.a $(TARGET)
	ranlib libnornir.a
	ar -M < buildlib.ar
install:
	mkdir -p $(NORNIR_PATH_INCLUDE)
	find . -name '*.hpp' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	cp ./libnornir.a $(NORNIR_PATH_LIB)/libnornir.a
uninstall:
	rm -r $(NORNIR_PATH_INCLUDE)/*
	rm $(NORNIR_PATH_LIB)/libnornir.a

