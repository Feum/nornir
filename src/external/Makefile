.PHONY = all

all:
# Checkout correct versions of external libraries
	mkdir -p mammut 
ifeq ("$(wildcard ./mammut/.git)","")
	git clone https://github.com/DanieleDeSensi/mammut.git #&& cd mammut && git checkout 47ef67d33a998279345ea6f8741b7fa9ec9847b7
endif
	mkdir -p knarr
ifeq ("$(wildcard ./knarr/.git)","")
	git clone https://github.com/DanieleDeSensi/knarr.git #&& cd knarr && git checkout aae8e80398a6a516705f0b5472db1e2881dae396
endif
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

