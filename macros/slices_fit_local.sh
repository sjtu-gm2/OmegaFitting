#!/bin/bash
# export LD_LIBRARY_PATH=$HOME/workspace/OmegaFitter/objs:${LD_LIBRARY_PATH}

version=oldfitter
clusterCut=Cut40 #Cut30 Cut40 Cut100
# nslice=${1}

inputROOT=/Users/cheng/workspace/Data/Run4U/${version}/hists_PU_${clusterCut}.root
inputLM=/Users/cheng/workspace/Data/Run4U/LM_Run4U.root
outputDir=./output/${version}_${clusterCut}

mkdir -p ${outputDir}

# echo $LD_LIBRARY_PATH
# 2 3 8 10


for nslice in 12
# for ((nslice=0;nslice<=19;nslice++))
do
	./main.exe ${inputROOT} ${inputLM} ${outputDir} ${nslice} | tee ./logs/${version}_${clusterCut}_${nslice}.log

	echo -e "FINISHED-${version}-${clusterCut}-${nslice}"
done

hadd -f ./output/${version}_${clusterCut}.root ./output/${version}_${clusterCut}/*.root

python ./EnergySlicesFitAna.py ${version} ${clusterCut}