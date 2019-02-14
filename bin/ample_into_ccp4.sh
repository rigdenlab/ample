#!/bin/bash
# Script for linking an ample checkout into the current CCP4 installation

ample_root=$( dirname "${BASH_SOURCE[0]}" )"/.."
ample_root=$(cd $ample_root && pwd -P)
ample_py=$ample_root/ample

ccp4_py=$CCP4/lib/py2/ample
ccp4_share=$CCP4/share/ample

BAK_EXT=".ccp4"

# Remove old symlinks and restore the original versions
if [ -L $ccp4_py ] ||  [ -L $ccp4_share ] ; then
    echo "Removing old symlinks and restoring backup"
	if [ -L $ccp4_py ]; then
	   rm  $ccp4_py
	   if [ -d ${ccp4_py}${BAK_EXT} ]; then
	     mv ${ccp4_py}${BAK_EXT} ${ccp4_py}
	   fi
	fi
	if [[ -L $ccp4_share ]]; then
	   rm  $ccp4_share
	   if [ -d ${ccp4_share}${BAK_EXT} ]; then
	     mv ${ccp4_share}${BAK_EXT} ${ccp4_share}
	   fi
	fi
	echo "Exiting after removing backup"
	exit
else
	echo "Linking in new AMPLE"
    mv  ${ccp4_py} ${ccp4_py}${BAK_EXT}
    mv ${ccp4_share} ${ccp4_share}${BAK_EXT}
    ln -s $ample_py $ccp4_py
    ln -s $ample_root $ccp4_share
    echo "Finished new links"
fi