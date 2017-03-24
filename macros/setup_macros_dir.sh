#!/bin/bash
make clean -f MakefileRIDF
make -f MakefileRIDF
make clean -f MakefileBeam
make -f MakefileBeam
make clean -f MakefileDCTPF
make -f MakefileDCTPF
make clean -f MakefileBDC
make -f MakefileBDC
make clean -f MakeRIDFoutput
make -f MakeRIDFoutput
make clean -f MakefileIC
make -f MakefileIC
