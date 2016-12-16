.PHONY = all

all:
	cd Mammut && git checkout 8d2c36cd20b653f7f58123b6b6ac5e5b69624e6b
	$(MAKE) -C Mammut
	$(MAKE) -C leo
	$(MAKE) -C queues msqueue
	$(MAKE) -C orlog
clean:
	$(MAKE) -C leo clean
	$(MAKE) -C queues clean
cleanall:
	$(MAKE) -C leo cleanall
	$(MAKE) -C queues cleanall
