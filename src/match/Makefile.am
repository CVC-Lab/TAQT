bindir = ../bin

bin_PROGRAMS = match contree readrawiv
match_SOURCES = main.cpp
match_LDADD = ../betti/libbetti.a ../par/libpar.a ../bio/libbio.a ../nr/libnr.a 

contree_SOURCES = contree.cpp
contree_LDADD = ../betti/libbetti.a ../par/libpar.a ../bio/libbio.a ../nr/libnr.a

readrawiv_SOURCES = readrawiv.cpp
readrawiv_LDADD = ../par/libpar.a ../bio/libbio.a ../nr/libnr.a 

INCLUDES = -I../par -I../betti -I../nr -I/home/junenim/Libs/XMLRPC/XmlRpcC -I..

