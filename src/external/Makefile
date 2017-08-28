.PHONY = all

all:
# Checkout correct versions of external libraries
	cd mammut && git checkout 4c609823f1ff268b1dd05210a93debbcc3bdfadb
	cd knarr && git checkout d418655ad7333d1886037cd94f76d3a31980990b
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
# Build graphite client
	$(MAKE) -C graphite-c-client all
clean:
	$(MAKE) -C leo clean
cleanall:
	$(MAKE) -C leo cleanall

