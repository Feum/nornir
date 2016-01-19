APPS                 = videoprocessing #blackscholes #qt-mandelbrot quicksort nqueens pbzip2

.PHONY: all clean cleanall

all:
#	cd qt-mandelbrot && qmake
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

