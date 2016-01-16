find_package ( Qt4 REQUIRED )
include( ${QT_USE_FILE} )

set ( CMAKE_CXX_FLAGS "--std=c++11")

set ( MANDSRCS
    main.cpp
    mandelbrotwidget.cpp
    renderthread.cpp
)

set ( MANDHDRS
        mandelbrotwidget.h
        renderthread.h
)

include_directories( 
    ${QT_INCLUDE_DIR}
    ${CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/atomic
    ../../src/external/Mammut/
)

QT4_WRAP_CPP( MOC_SRCS ${MANDHDRS} )
 
find_package (Threads)
add_executable ( mandelbrot_ff ${MANDSRCS} ${MOC_SRCS} )
target_link_libraries ( mandelbrot_ff ${QT_LIBRARIES} )
target_link_libraries (mandelbrot_ff ${CMAKE_THREAD_LIBS_INIT})
