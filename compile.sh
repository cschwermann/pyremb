#!/bin/bash

gfortran -DLIBXC -o Test Precision_mod.F90 Output_mod.F90 Types_mod.F90 System_mod.F90 Utils_mod.F90 external/Xcpot_libxc_mod.F90 Xcpot_mod.F90 Test.F90 -lxc -I/usr/include/ -I./external

