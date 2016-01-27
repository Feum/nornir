# Makefile
#

# Define MANDEL_NO_DISPLAY if you do not want to display image
# ex: export MANDEL_NO_DISPLAY=1  or
#     setenv MANDEL_NO_DISPLAY 1
#

MANDEL_NO_DISPLAY=1

ifdef MANDEL_NO_DISPLAY
DIS=-DNO_DISPLAY
else
DIS=
endif


CXXFLAGS             += -DNO_CMAKE_CONFIG
CFLAGS               += $(DIS) 
LDFLAGS              += $(DIS) -L/usr/X11R6/lib -L.
INCS                 += -I. -I/usr/X11R6/include
LDLIBS                 += -lpthread -lmarX2 -lX11 -lm

TARGET               = mandel_ff mandel_seq 
SOURCES		     = mandel_ff.cpp mandel_seq.c marX2.c
SRC = .
INC = .

.PHONY: all clean cleanall

all:$(TARGET)

libmarX2.a: $(SRC)/marX2.c $(INC)/marX2.h
	gcc -c $(OPTIMIZE_FLAGS) $(CFLAGS) -I/usr/X11R6/include $(INCS) $(SRC)/marX2.c 
	ar -rv libmarX2.a marX2.o
	ranlib libmarX2.a

mandel_seq: mandel_seq.c libmarX2.a
	$(CC)  $(CFLAGS) $(OPTIMIZE_FLAGS) $(CFLAGS)  $(INCS) $(LDFLAGS)  -o mandel_seq mandel_seq.c $(LDLIBS)

mandel_ff: mandel_ff.cpp  libmarX2.a
	$(CXX)  $(CXXFLAGS) $(OPTIMIZE_FLAGS)  $(CFLAGS) $(INCS) $(LDFLAGS)  -o mandel_ff mandel_ff.cpp $(LDLIBS)

clean:
	rm -f libmarX2.a *.o

cleanall: clean
	rm -fr $(TARGET) *~ *.dSYM


