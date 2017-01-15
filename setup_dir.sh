g#!/bin/bash

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
PS3="Please choose your working environment:"
select option in S015 Fishtank RIKENHPC MSUHPC other
do
    case $option in
        S015)
            echo "S015";;
        Fishtank)
            echo "Fishtank";;
        RIKENHPC)
            echo "RIKENHPC";;
        MSUHPC)
            echo "MSUHPC";;
        other)
            echo "Please link ridf folder manually";;
     esac
done
