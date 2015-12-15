TARGET = knob.o node.o parameters.o predictors.o configuration.o utils.o

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
	ar rc libgrape.a $(TARGET)
	ranlib libgrape.a
	ar -M < buildlib.ar
install:
	mkdir -p $(GRAPE_PATH_INCLUDE)
	find . -name '*.hpp' -exec cp --parents \{\} $(GRAPE_PATH_INCLUDE) \;
	cp ./libgrape.a $(GRAPE_PATH_LIB)/libgrape.a
uninstall:
	rm -r $(GRAPE_PATH_INCLUDE)/*
	rm $(GRAPE_PATH_LIB)/libgrape.a
