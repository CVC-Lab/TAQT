CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (CVC_HelloWorld)

# Global Variables

#SET (CMAKE_CXX_FLAGS "-Wno-deprecated -Wno-non-template-friend")
if(CMAKE_HOST_UNIX)
  SET (CMAKE_CXX_FLAGS $ENV{CPPFLAGS})
  SET (CMAKE_LD_FLAGS $ENV{LDFLAGS})
endif(CMAKE_HOST_UNIX)

#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

ADD_SUBDIRECTORY (src)
