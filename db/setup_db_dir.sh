#!/bin/bash
#this script will create links to online analysis configuration files. These links can be changed by other bash scripts
if [ ! -e BigRIPSIC.xml ]; then
  ln -sf ./BigRIPSIC/BigRIPSIC.xml BigRIPSIC.xml
fi
if [ ! -e BigRIPSPlastic.xml ]; then
ln -sf ./BigRIPSPlastic/BigRIPSPlastic.xml BigRIPSPlastic.xml
fi
if [ ! -e BigRIPSPPAC.xml ]; then
ln -sf ./BigRIPSPPAC/BigRIPSPPAC.xml BigRIPSPPAC.xml
fi
if [ ! -e NEBULA.xml ]; then
ln -sf ./NEBULA/NEBULA.xml NEBULA.xml
fi
if [ ! -e NEULAND.xml ]; then
ln -sf ./NEULAND/NEULAND.xml NEULAND.xml
fi
if [ ! -e NEULANDVETO.xml ]; then
ln -sf ./NEULANDVETO/NEULANDVETO.xml NEULANDVETO.xml
fi
if [ ! -e SAMURAIBDC1.xml ]; then
ln -sf ./SAMURAIBDC/SAMURAIBDC1.xml SAMURAIBDC1.xml
fi
if [ ! -e SAMURAIBDC2.xml ]; then
ln -sf ./SAMURAIBDC/SAMURAIBDC2.xml SAMURAIBDC2.xml
fi
if [ ! -e SAMURAIBPC.xml ]; then
ln -sf ./SAMURAIBPC/SAMURAIBPC_nothing.xml SAMURAIBPC.xml
fi
if [ ! -e SAMURAIFDC1.xml ]; then
ln -sf ./SAMURAIFDC/SAMURAIFDC1.xml SAMURAIFDC1.xml
fi
if [ ! -e SAMURAIFDC2.xml ]; then
ln -sf ./SAMURAIFDC/SAMURAIFDC2.xml SAMURAIFDC2.xml
fi
if [ ! -e SAMURAIHOD.xml ]; then
ln -sf ./SAMURAIHOD/SAMURAIHOD.xml SAMURAIHOD.xml
fi
if [ ! -e SAMURAIPlastic.xml ]; then
ln -sf ./SAMURAIPlastic/SAMURAIPlastic.xml SAMURAIPlastic.xml
fi
