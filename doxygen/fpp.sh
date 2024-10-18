#!/bin/bash
#cpp -lang-fortran -traditional-cpp -D_LANGUAGE_FORTRAN -I/mnt/appl/x86_64/petsc/3.3-p7/include/ -I/mnt/appl/x86_64/mpich/mpich-3.0.4_ifc2013/include -I/mnt/appl/x86_64/slepc/3.3-p4/include/ -I/mnt/appl/x86_64/slepc/3.3-p4/include/finclude/ $* | grep -v -e '#' | grep -v '^$'
#cat $*
ifort -fpp -E -I/mnt/appl/x86_64/petsc/3.12.0/include/ -I/mnt/appl/x86_64/mpich/mpich-3.3.1/include -I/mnt/appl/x86_64/slepc/3.12.0/include/ $* | grep -v -e '#' | grep -v '^$'
