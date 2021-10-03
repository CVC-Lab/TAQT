CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/match/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/nr
  ../../inc/par
  ../../inc/betti
  ../../inc/XmlRpcC
  ../../inc/bio
)

ADD_EXECUTABLE (match main.cpp)
TARGET_LINK_LIBRARIES(match
	betti
	par
	bio
	nr
	XmlRpcC
)


ADD_EXECUTABLE (readrawiv readrawiv.cpp)

TARGET_LINK_LIBRARIES(readrawiv
  par
  bio
)


ADD_EXECUTABLE (contree contree.cpp)

TARGET_LINK_LIBRARIES(contree
  betti
  par
  bio
  nr
)
