#!/bin/bash


export OMP_NUM_THREADS=4
prmfile=2dir.prm

./2dir.exe $prmfile > 2dir.out 2>&1
