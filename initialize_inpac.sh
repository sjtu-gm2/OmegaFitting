# setup root and c++

main=Main

if [[ ${1} -eq "run23" ]];then
    run=2.0
elif [[ ${1} -eq "run4" ]];then
    run=4.0
else
    echo "unknown dataperiod ${1}" 
    exit 1
fi

echo -e "initialize framework for Run-${run} analysis"

source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

mkdir ./run && mkdir ./build

cd OmegaFitting && ln -s $PWD/launch.py $PWD/example.json ../run

cd ../build && cmake ../OmegaFitting -DRun=${run} -DMain=${main} -DUSE_JSON=OFF && cmake --build . && cd ../run && mkdir logs


