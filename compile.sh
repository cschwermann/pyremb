#!/bin/bash

rm *.mod Test Do_test

gfortran -DLIBXC -g -fbacktrace -ffpe-summary=denormal,zero,invalid,overflow,underflow -o Test Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 Peremb_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Test.F90 -lxc -I/usr/include/ -I./external

f90wrap -m peremb Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Peremb_mod.F90 

f2py-f90wrap -c -m _peremb OBJ_FILES f90wrap_*.f90 *.o -lxc -I/usr/include/ -I./external -I. -DLIBXC


#nagfor -o Test Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Test.F90 -lxc -I/usr/include/ -I./external

#documentation
rm -r docs/*
ford Pyremb-project.md
