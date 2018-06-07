#!/bin/sh
#This is the script that is supposed to start a parallel run on a machine
#with two quad core processors. This is just from command line, so no 
#scheduling (e.g. PBS) is used.

#Setting enviromenemt variable from Intel fortran compiler and openmpi
module add intel/compiler/64/11.1/046 openmpi/intel/64/1.6.3 netcdf/intel/64/4.1.3

export LD_LIBRARY_PATH=/opt/netcdf_intel/4.1.3/64/lib:/cm/shared/apps/gcc/4.7.0/lib:/cm/shared/apps/gcc/4.7.0/lib64:/cm/shared/apps/torque/3.0.5/lib/:/cm/shared/apps/maui/3.3.1/lib

export LD_LIBRARY_PATH=/opt/python/gcc/64/2.7.3/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/netcdf_intel/4.1.3/64/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/openmpi/intel/64/1.6.3/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cm/shared/apps/intel/Compiler/11.1/046/lib/intel64/:$LD_LIBRARY_PATH

#The name of the telemac cas file
telmod=$1
casfile=$2

#Specify telemac root folder and the corresponding folders for teh configs,
#builds and python scripts, respectively
TELEMAC=/opt/telemac/v6p3r1
SYSTELCFG=$TELEMAC/configs
TELBUILD=$TELEMAC/builds/ifortmpicdf
TELPYT=$TELEMAC/scripts/python27

#name of the config file and the option to be used (ifortmpi), respectively 
cfgfile=systel_WL.cfg
cfopt=ifortmpicdf

export PATH=/opt/python/gcc/64/2.7.3/bin:$TELBUILD/bin:$TELPYT:$PATH
export LD_LIBRARY_PATH=$TELBUILD/lib:/opt/metis/intel/64/5.1.0/lib/:$LD_LIBRARY_PATH

echo " "
#This is the actual call to run telemac2d
python $TELPYT/runcode.py $telmod -c $cfopt -f $SYSTELCFG/$cfgfile -s $casfile
