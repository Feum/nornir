TARGET = manager.o knob.o node.o parameters.o predictors.o utils.o

.PHONY: all clean cleanall install uninstall lib

all:
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
install:
	mkdir -p $(GRAPE_PATH_INCLUDE)
	find . -name '*.hpp' -exec cp --parents \{\} $(GRAPE_PATH_INCLUDE) \;
	cp ./libgrape.a $(GRAPE_PATH_LIB)/libgrape.a
uninstall:
	rm -r $(GRAPE_PATH_INCLUDE)/*
	rm $(GRAPE_PATH_LIB)/libgrape.a

