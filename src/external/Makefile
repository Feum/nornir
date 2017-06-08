.PHONY = all

all:
# Build mammut
	$(MAKE) -C mammut
# Build leo
ifeq ("$(wildcard ./leo/leo.o)","")
	$(MAKE) -C leo
endif
# Build queues
ifeq ("$(wildcard ./queues/hzdptr.o)","")
	$(MAKE) -C queues msqueue
endif
# Build knarr
	$(MAKE) -C knarr

clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

