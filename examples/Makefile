MAMMUTROOT           = /home/desensi/Mammut
export INCS          = -I$(realpath ..) -I$(MAMMUTROOT)
export LIBS          = -lpthread -L$(MAMMUTROOT)/mammut -lmammut -lprotobuf-lite -lrt
APPS                 = cholesky pbzip2 matmul

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

