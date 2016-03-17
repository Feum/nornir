APPS                 = canneal simple_mandelbrot blackscholes pbzip2 videoprocessing quicksort 

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

