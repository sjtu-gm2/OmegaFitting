#!/bin/bash
# source ./OmegaFitting/compile.sh -t ./OmegaFitting/Run4U_nominal.cpp

# export LD_LIBRARY_PATH=$PWD/objs:${LD_LIBRARY_PATH}

version=Run4U_Tmethod_Double

EXE=Run4U_nominal
startTimeShift=0
seed=0
chainFit=1

inputROOT=/Users/cheng/WorkRun4/data/nominal_fit/Run4U_empirical_nominal.root
inputLM=/Users/cheng/WorkRun4/data/nominal_fit/Run4U_lm_int_sjtu.root

wiggle=1 # is triple?

outputDir=./output/${version}




mkdir -p ${outputDir}
mkdir -p ./logs/${version}
mkdir -p ./json/${version}

for ((startTimeShift=2;startTimeShift<=2;startTimeShift++)) 
do
    tag=Seed${seed}_Time${startTimeShift}
    InitValues=./run4_baseline_T.json
    # InitValues=/Users/cheng/WorkRun4/json/Run4U_Tmethod_nominal/Seed0_Time0.json
    ./build/MAIN ${inputROOT} ${wiggle} ${inputLM} ${outputDir} ${InitValues} ${startTimeShift} ${version}_${tag} ${chainFit} | tee ./logs/${version}/${tag}.log
done
