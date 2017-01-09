#!/bin/bash

if [ -e $TARTSYS ]; then
  if [ -e src ]; then
    echo "src exists already"
  fi
  if [ ! -e src ]; then
    echo "using $TARTSYS as src directory"
    ln -sf $TARTSYS src
  fi
else
echo "cannot find TARTSYS location. Check that anaroot is installed and sourced."
fi
echo "attempting to setup initial db links"
cd db/
./setup_db_dir.sh
