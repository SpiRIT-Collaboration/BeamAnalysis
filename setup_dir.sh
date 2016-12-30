#!/bin/bash

if [ -e $TARTSYS ]; then
  if [ -e src ]; then
    echo "src exists already"
  fi
  if [ ! -e src ]; then
    echo "using $TARTSYS as src directory"
    ln -sf $TARTSYS src
  fi
fi
