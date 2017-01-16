#!/bin/bash

#two inputs are used, the first and the last run that we want to run RIDFtoROOT for
#We also need the input file, ridf_events.csv
START=$1
END=$2

#check for input runs and adjust accordingly

if [ $# -eq 0 ]; then
echo "no run number provided, please provide the first run you want to use:"
read inRun
echo "please provide the last run you want to use:"
read inRun2
START="$inRun"
END="$inRun2"
fi

if [ $# -eq 1 ]; then
END="$START"
echo "evaluating only run "$END
fi


#read through ridf_events.csv and allocate each column
while IFS=$'\n', read -r -a input;
do
RUN=${input[0]}
EVENTS=${input[1]}
RUNTYPE=${input[2]}
IC=${input[3]}
PPAC=${input[4]}
PLAS=${input[5]}
NEU=${input[6]}
NEUVETO=${input[7]}
BDC1=${input[8]}
BDC2=${input[9]}
TOF=${input[10]}
DCTPF=${input[11]}

if  [ "$RUN" -eq "$RUN" ] 2>/dev/null; then
if [ "$RUN" -ge $START ] && [ "$RUN" -le $END ] && [ "$EVENTS" -ge 1 ]; then
echo "$RUN"
cd db/
./setup_db_dir.sh
cd ../
if [ ${#IC} -ge 1 ]; then
ln -sf BigRIPSIC/BigRIPSIC."$IC".xml ../db/BigRIPSIC.xml
else
ln -sf BigRIPSIC/BigRIPSIC.xml ../db/BigRIPSIC.xml
fi
if [ ${#PLAS} -ge 1 ]; then
ln -sf BigRIPSPlastic/BigRIPSPlastic."$PLAS".xml ../db/BigRIPSPlastic.xml
else
ln -sf BigRIPSPlastic/BigRIPSPlastic.xml ../db/BigRIPSPlastic.xml
fi
if [ ${#PPAC} -ge 1 ]; then
ln -sf BigRIPSPPAC/BigRIPSPPAC."$PPAC".xml ../db/BigRIPSPPAC.xml
else
ln -sf BigRIPSPPAC/BigRIPSPPAC.xml ../db/BigRIPSPPAC.xml
fi
if [ ! -e dctpf/dc_tpf_$DCTPF.root ]; then
  echo "dc_tpf file not found, creating it"
  ./SAMURAIDCTPF $DCTPF
fi
rm dctpf/dc_tpf.root
ln -sf dctpf/dc_tpf_$DCTPF.root dctpf/dc_tpf.root
if [ ${#TOF} -ge 1 ]; then
./RIDFtoROOT $RUN $TOF
else
echo "no TOF information in ridf_events.csv for run "$RUN", using 280 ns"
./RIDFtoROOT $RUN 280.
fi

else
:
fi
#end filter for run number

else
echo "bypassing detected header"
fi
#end filter for header

done < "ridf_events.csv"
