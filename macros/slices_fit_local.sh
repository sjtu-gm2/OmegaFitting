#!/bin/bash
# export LD_LIBRARY_PATH=$HOME/workspace/OmegaFitter/objs:${LD_LIBRARY_PATH}

version=newfitter
clusterCut=All #Cut30 Cut40 Cut100
EXE=energy_slices_scan
InitValues=initial_values_${version}_${clusterCut}.json

inputROOT=/Users/cheng/workspace/Data/Run4U/${version}/hists_PU_${clusterCut}.root
inputLM=/Users/cheng/workspace/Data/Run4U/LM_Run4U.root
outputDir=./output/${version}_${clusterCut}

mkdir -p ${outputDir}

for nslice in 3 13 19
# for ((nslice=0;nslice<=19;nslice++))
do
	./${EXE}.exe ${inputROOT} ${inputLM} ${outputDir} ${nslice} ${InitValues} | tee ./logs/${version}_${clusterCut}_${nslice}.log

	echo -e "FINISHED-${version}-${clusterCut}-${nslice}"
done

hadd -f ./output/${version}_${clusterCut}.root ./output/${version}_${clusterCut}/*.root

python ./srcs/macros/EnergySlicesFitAna.py ${version} ${clusterCut}
python ./srcs/macros/ExtractParameters.py ${version} ${clusterCut}