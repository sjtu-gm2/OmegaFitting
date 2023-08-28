# setup root and c++

main=Main
first_build=${2:-true}

if [[ ${1} == "run23" ]];then
    run=2.0
elif [[ ${1} == "run4" ]];then
    run=4.0
else
    echo "unknown dataperiod ${1}" 
    exit 1
fi

source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

if [[ ${first_build} == true ]]
then
    echo -e "initialize framework for Run-${run} analysis"

    mkdir ./run && mkdir ./build
    cd OmegaFitting && ln -s $PWD/launch.py $PWD/example.json $PWD/configs ../run
    cd ../build && cmake ../OmegaFitting -DRun=${run} -DMain=${main} -DUSE_JSON=OFF && cmake --build . && cd ../run && mkdir logs
else
    echo -e "rebuild framework for Run-${run} analysis"
    
    mkdir -p build
    cd build && cmake ../OmegaFitting -DRun=${run} -DMain=${main} -DUSE_JSON=OFF && cmake --build . && cd ..
fi



