#!/bin/bash

version=${1}
clusterCut=${2}
nslice=${3}

source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
export LD_LIBRARY_PATH=/home/chencheng/OmegaFitter/objs:${LD_LIBRARY_PATH}

inputROOT=/home/chencheng/data/Run4U/${version}/hists_PU_${clusterCut}.root
inputLM=/home/chencheng/data/Run4U/LM_Run4U.root
outputDir=/home/chencheng/OmegaFitter/output/${version}_${clusterCut}



./main.exe ${inputROOT} ${inputLM} ${outputDir} ${nslice}

echo -e "FINISHED-${version}-${clusterCut}-${nslice}"
