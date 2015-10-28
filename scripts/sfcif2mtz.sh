#!/bin/bash

inCif=$1
name=`echo $1 | cut -c -4`

inPdb=$2
cell=`cat $2 | grep CRYST1 | cut -c 7-54`
echo $cell
symm=`cat $2 | grep CRYST1 | cut -c 56-66 | sed 's/ //g'`
echo $symm

cif2mtz hklin $inCif hklout ${name}_sf.mtz << eof
cell $cell
symm $symm
end
eof

#ctruncate -mtzin "${name}_sf.mtz" -mtzout "${name}_ctruncate.mtz" -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" -Imean
ctruncate -mtzin "${name}_sf.mtz" -mtzout "${name}_ctruncate.mtz" -colin "/*/*/[I,SIGI]"

uniqueify "${name}_ctruncate.mtz"

