CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES 
  actree.cpp contourtree.cpp interval_tree.cpp
  arc.cpp disjointset.cpp node.cpp volume3dcritical.cpp
  Volume2DCritical.cpp dual.cpp multicon.cpp moment.cpp
  tetinteg.cpp VolCritical.cpp
)

FILE(GLOB INCLUDE_FILES ../../inc/betti/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/betti
  ../../inc/par
  ../../inc/nr
)

ADD_LIBRARY (betti STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
