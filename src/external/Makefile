.PHONY = all

all:
ifneq ("$(wildcard Mammut/mammut/libmammut.a)","")
	$(MAKE) -C Mammut
endif
ifneq ("$(wildcard leo/leo.o)","")
	$(MAKE) -C leo
endif
ifneq ("$(wildcard queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
ifneq ("$(wildcard ./orlog/src/liborlog.a)","")
	$(MAKE) -C orlog
endif
clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

