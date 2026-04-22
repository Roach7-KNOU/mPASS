# Library versioning
set(SOVERSION 4)
set(LIBRARY_VERSION 4.0.1)
project (pass)

set(SRCS
    pass.cpp
)
set(LIBSRCS
    reader.cpp cal.cpp soga.cpp
)

add_library(libs SHARED ${LIBSRCS})
add_executable(pass ${SRCS})

target_link_libraries ( pass libs)

ADD_DEFINITIONS( -DSQLITE_THREADSAFE=1 -g -O2 )
