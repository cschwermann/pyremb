#!/bin/bash

gfortran -DLIBXC -g -fbacktrace -ffpe-summary=denormal,zero,invalid,overflow,underflow -o Test Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Test.F90 -lxc -I/usr/include/ -I./external


#nagfor -o Test Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Test.F90 -lxc -I/usr/include/ -I./external
