#!/bin/bash

if [ -e src ]; then
echo "src exists"
fi
if [ ! -e src ]; then
echo "Using $SIMPATH as src directory"
ln -sf $SIMPATH src
fi

