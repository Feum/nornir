TARGET  = blackscholes
M4_BASE = .
MACROS  = c.m4.pthreads

ifdef source
	ifeq "$(source)" "simd" 
		SRC = blackscholes.simd.c
		CXXFLAGS += -msse3
	endif
else
	SRC	= blackscholes.c
endif

CSRC    = $(SRC)

version = ff
#chk_err = 1

ifdef version
  ifeq "$(version)" "pthreads"
    M4_SRC    = blackscholes.m4.cpp
    CSRC      = $(M4_SRC)
    MT        = -DENABLE_THREADS
    CXXFLAGS += -pthread
  endif
  ifeq "$(version)" "openmp"
    MT        = -DENABLE_OPENMP
  endif
  ifeq "$(version)" "tbb"
    MT        = -DENABLE_TBB
  endif
  ifeq "$(version)" "ff"
    CXXFLAGS += --std=c++11 -O3
    MT        = -DENABLE_FF
    LDLIBS      += -pthread
  endif
else
  MT        = 
endif



# Default build single precision version
NCO     = -DNCO=4

ifdef chk_err
ERR     = -DERR_CHK
endif

ifdef single
NCO = -DNCO=4
endif

ifdef size
SZ = -DN=$(size)
else
SZ = -DN=960
endif

ifdef double
NCO = -DNCO=2
endif

CXXFLAGS += $(MT) $(NCO) $(FUNC) $(ERR)

all: $(TARGET)

$(TARGET): $(CSRC)
	$(CXX) $(CXXFLAGS) $(CSRC) $(INCS) $(LDFLAGS) $(LDLIBS) -o $(TARGET)
	$(CXX) $(CXXFLAGS) inputgen.c $(INCS) $(LDFLAGS) $(LDLIBS)  -o inputgen

$(M4_SRC): $(SRC)
	$(M4) $(M4_BASE)/$(MACROS) $< > $@

clean:
	rm -f $(TARGET) $(M4_SRC)
cleanall:
	make clean
