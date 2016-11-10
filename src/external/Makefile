.PHONY = all

all:
ifeq ("$(wildcard Mammut/mammut/libmammut.a)","")
	$(MAKE) -C Mammut
endif
ifeq ("$(wildcard leo/leo.o)","")
	$(MAKE) -C leo
endif
ifeq ("$(wildcard queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
ifeq ("$(wildcard ./orlog/src/liborlog.a)","")
	$(MAKE) -C orlog
endif
clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

