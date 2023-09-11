# setup root and c++

main=Main
first_build=${1:-true}

source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

if [[ ${first_build} == true ]]
then
    echo -e "initialize framework for omega_a analysis"

    mkdir ./run && mkdir ./build
    cd OmegaFitting && ln -s $PWD/launch.py $PWD/example.json $PWD/configs ../run
    cd ../build && cmake ../OmegaFitting -DMain=${main} -DUSE_JSON=OFF && cmake --build . && cd ../run && mkdir logs
else
    echo -e "rebuild framework for omega_a analysis"
    
    mkdir -p build
    cd build && cmake ../OmegaFitting -DMain=${main} -DUSE_JSON=OFF && cmake --build . && cd ..
fi



