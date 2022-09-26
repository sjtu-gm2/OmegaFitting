# setup root and c++
run=2.0
main=Run2

echo -e "initialize framework for Run-${run} analysis"


source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

mkdir ./run && mkdir ./build
cd run && ln -s launch_fitter_run2.py run_scan.sh *.condor ../run

cd ../build && cmake ../OmegaFitting -DRun=${run} -DMain=Run2 -DUSE_JSON=OFF && cmake --build .


