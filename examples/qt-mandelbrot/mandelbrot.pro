INCLUDEPATH +=  ../../src/external/fastflow/
INCLUDEPATH +=  ../../src/external/Mammut/
INCLUDEPATH +=  /usr/include/libxml
VPATH +=  ../../src/external/fastflow/
LIBS +=  -L ../../src -lgrape -pthread -lrt -lm -lmlpack -llapack -lblas -lgsl -lgslcblas
QMAKE_CXXFLAGS += -std=c++0x
SOURCES += main.cpp
SOURCES += mandelbrotwidget.cpp
SOURCES += renderthread.cpp  
HEADERS += renderthread.h
HEADERS += mandelbrotwidget.h
HEADERS += renderthread.h
HEADERS += ff/farm.hpp
HEADERS += ff/node.hpp

TARGET = mandelbrot-qt


