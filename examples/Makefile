APPS                 = simple-mandelbrot blackscholes pbzip2 videoprocessing #blackscholes quicksort nqueens pbzip2

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

