MAMMUTROOT           = /home/daniele/Mammut
export INCS          = -I$(realpath ..) -I$(MAMMUTROOT) -I/usr/include/libxml2
export LIBS          = -lpthread -L$(MAMMUTROOT)/mammut -lmammut -lprotobuf-lite -lrt -lmlpack -llapack -lblas
APPS                 = pbzip2  #cholesky matmul

.PHONY: all clean cleanall

all:
	for APP in $(APPS) ; do \
		$(MAKE) -C $$APP ; \
	done
clean: 
	for APP in $(APPS) ; do \
		$(MAKE) -C $$APP clean; \
	done
cleanall:
	for APP in $(APPS) ; do \
		$(MAKE) -C $$APP cleanall; \
	done

