#!/bin/bash

if [ -e $TARTSYS ]; then
  if [ -e src ]; then
    echo "src link exists already"
  fi
  if [ ! -e src ]; then
    echo "Using $TARTSYS as src directory"
    ln -sf $TARTSYS src
  fi
else
echo "Cannot find TARTSYS location. Check that anaroot is installed and sourced."
fi
echo "Creating initial db links"
cd db/
./setup_db_dir.sh
echo "please create ridf link to folder containing ridf files:"
echo "ln -sf RIDF/FILE/LOCATION ridf"
cd ../macros
./setup_macros_dir.sh
cd ../
if [ ! -e ridf ]; then
PS3="Please choose your working environment:"
select option in S015 Fishtank RIKENHPC MSUHPC other
do
    case $option in
        S015)
            ln -sf /home/s015/ridf/sdaq02 ridf
            break;;
        Fishtank)
            ln -sf /mnt/spirit/rawdata/ridf ridf
            break;;
        RIKENHPC)
            ln -sf /data/Q16264/rawdata/ridf/sdaq02 ridf
            break;;
        MSUHPC)
            ln -sf /mnt/research/spirit/SPIRIT_TPC/data/ridf ridf
            break;;
        other)
            echo "Please link ridf folder manually"
            break;;
     esac
done
else
    echo "ridf link exists, leaving unchanged"
fi
