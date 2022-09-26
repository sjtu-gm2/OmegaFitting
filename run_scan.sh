#! /bin/bash
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

method=${1}
random_type=${2}
seed_num=${3}
scan_type=${4}

cd /home/chencheng/Fitter_Wa/condor
python ./launch_fitter_condor.py ${method} ${random_type} ${seed_num} ${scan_type}