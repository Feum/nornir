.PHONY = all

all:
# Build mammut
ifeq ("$(wildcard ./Mammut/mammut/libmammut.a)","")
	$(MAKE) -C Mammut
endif
# Build leo
ifeq ("$(wildcard ./leo/leo.o)","")
	$(MAKE) -C leo
endif
# Build queues
ifeq ("$(wildcard ./queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
# Build orlog
ifeq ("$(wildcard ./orlog/src/liborlog.a)","")
	$(MAKE) -C orlog
endif

clean:
	$(MAKE) -C leo clean
	$(MAKE) -C queues clean
cleanall:
	$(MAKE) -C leo cleanall
	$(MAKE) -C queues cleanall
