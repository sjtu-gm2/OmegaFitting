#! /bin/bash
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
timeShift=${1}

cd /home/chencheng/Fitter_Wa/condor
python ./launch_fitter_condor.py ${timeShift}