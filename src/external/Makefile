.PHONY = all

all:
# Build mammut
ifeq ("$(wildcard ./mammut/mammut/libmammut.a)","")
	$(MAKE) -C mammut
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
cleanall:
	$(MAKE) -C leo cleanall

