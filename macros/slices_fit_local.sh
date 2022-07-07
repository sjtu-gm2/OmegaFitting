#!/bin/bash
export LD_LIBRARY_PATH=$HOME/workspace/OmegaFitter/objs:${LD_LIBRARY_PATH}

version=Run4U
clusterCut=OldReco
EXE=EScan
endTime=350

InitValues=fitted_values_Run4URatioCutoffScan_0.json
# InitValues=initial_values_newfitter_All.json

inputROOT=/Users/cheng/workspace/Data/Run4U/RatioCutoff_40MeV/ratioCutOff_v10_03_01/PUCorrected_Run4U_v10_03_01_ratioCutold_gridRun.root
inputLM=/Users/cheng/workspace/Data/Run4U/LM_Run4U.root
outputDir=./output/${version}_${clusterCut}

mkdir -p ${outputDir}

# for nslice in 8 16
for ((nslice=0;nslice<=19;nslice++))
do
	./${EXE}.exe ${inputROOT} ${inputLM} ${outputDir} ${nslice} ${InitValues} ${endTime}| tee ./logs/${version}_${clusterCut}_${nslice}.log

	echo -e "FINISHED-${version}-${clusterCut}-${nslice}"
done

hadd -f ./output/${version}_${clusterCut}.root ./output/${version}_${clusterCut}/*.root

python ./srcs/macros/EnergySlicesFitAna.py ${version} ${clusterCut}
python ./srcs/macros/ExtractParameters.py ${version} ${clusterCut}

# version=newfitter
# clusterCut=Cut40 #Cut30 Cut40 Cut100

# python ./srcs/macros/EnergySlicesFitAna.py ${version} ${clusterCut}

# 0 1 3 5 10 11 12 13 14 15 16 19 
# 0 3 10 12 13 15 16 19
# 19 16 15 13
# 14
# 11 12 16 18
# 18 16 15 12 10 5 3 0
# 6 5 3 0
# 0 1 2 4 5 6 8 10 11 12 14 16 17 19
# 1 4 5 6 7 12 13 14 16 17 18 19
# 6 14 15 16 17 19
# 0 1 2 3 4 5 7 8 9 10 11 12 13 18
# 2 3 13 19
# 0 1 4 5 6 7 8 9 10 11 12 14 15 16 17 18