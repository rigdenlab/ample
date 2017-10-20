#!/bin/bash

#
# Script to take a ample repository and insert it into the current CCP4 installation
#

if [ -z $CCP4 ]; then
    echo "CCP4 variable is not set - cannot do anything!"
    exit
fi

ample_topdir=`(cd ..; pwd)`
ample_srcdir=`(cd ../ample; pwd)`

ccp4_topdir=$CCP4/share/ample
ccp4_srcdir=$CCP4/lib/py2/ample

if [ ! -d $ccp4_topdir ] || [ ! -d $ccp4_srcdir ]; then
    echo "Cannot find CCP4 directories $ccp4_topdir and $ccp4_srcdir !"
    exit
fi

echo "** Symlinking into CCP4 directories $ccp4_topdir and $ccp4_srcdir **"
if [ -L $ccp4_topdir ]; then
    rm $ccp4_topdir
elif [ -d $ccp4_topdir ] && [ ! -d  ${ccp4_topdir}.bak ]; then
    mv $ccp4_topdir ${ccp4_topdir}.bak
fi
ln -s $ample_topdir $ccp4_topdir


if [ -L $ccp4_srcdir ]; then
    rm $ccp4_srcdir
elif [ -d $ccp4_srcdir ] && [ ! -d  ${ccp4_srcdir}.bak ]; then
    mv $ccp4_srcdir ${ccp4_srcdir}.bak
fi
ln -s $ample_srcdir $ccp4_srcdir

# GUI2
gui2dir=$CCP4/share/ccp4i2/wrappers/AMPLE
if [ ! -d $gui2dir ]; then
	mkdir -p $gui2dir
else
	mv $gui2dir/script ${gui2dir}/script.bak
fi
ln -s $ample_topdir/i2 $gui2dir/script
cp $ample_topdir/i2/AMPLE-icon.svg   $CCP4/share/ccp4i2/qticons/AMPLE.svg
