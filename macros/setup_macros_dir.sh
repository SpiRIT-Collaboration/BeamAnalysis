#!/bin/bash
make clean -f MakefileRIDF
make -f MakefileRIDF
make clean -f MakefileBeam
make -f MakefileBeam
make clean -f MakefileDCTPF
make -f MakefileDCTPF
