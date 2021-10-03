#!/bin/sh

#
# build.sh
# arand, 4-25-2011
#
# Script for the Home NMI Build and Test Lab (BaTLab)
# https://nmi.cs.wisc.edu/
#
# This runs the build process and saves the resulting binaries
# in a tar directory.
#

cmake ./
make

tar -c -z -f ./results.tar.gz ./bin 
# important: the file must be named results.tar.gz
#            as that is the only file which is saved
#            by the BaTlab software