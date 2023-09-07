# Omega-a fitter for g-2 analysis

## Introduction
Perform omega-a fitting on the wiggle histogram. This fitter have features:

1. Integrated Omega-a `Blinder`
2. Easily submit fitting jobs to `inpac condor` for seed scan, start time scan, adhoc scan, etc.


To use this fitter, prepare the ROOT files with

1. Wiggle histogram 
1. Lost muon integral histogram
1. Initial values saved as `TArrayD` array or `TF1` function

Then write up the configuration `.json` file with necessary fitting setup and submit the fitter to `condor (inpac)` or run `local`


## An example with inpac-condor

Clone git and build with Run4 `blinding`

```
mkdir Fitter_wa && cd Fitter_wa
git clone git@github.com:sjtu-gm2/OmegaFitting.git
source ./OmegaFitting/initialize_inpac.sh run4
```

The `initialize_inpac.sh` script will build the project as well as create a `run` directory and link some configuration files into it.

If you already have those files and just want to rebuild the project, you can run
```
source ./OmegaFitting/initialize_inpac.sh run4 false
```
in the `Fitter_wa` directory.

**$\rm\color{MediumPurple}{[Note]}$** Whenever edtting the source file in the `Fitter` direectory, you have to build again by
```
cmake --build ./build
```
in the `Fitter_wa` directory.

Now, you can find an example of configuration file, `example.json`. You can perform several tasks with this configuration file:

1. Energy bin scan:
```
# local run
python ./launch.py --run example.json energy_bin dummy Run4A 0
# submit to condor
python ./launch.py --submit example.json energy_bin
```
Lines starting with `#` are comments. The option `--submit` will luanch `ceil(525 / 10) = 53` jobs in total to `inpac condor`.

2. Start time scan
```
# local run
python ./launch.py --run example.json start_time_scan T_1700 Run4A 0
# submit to condor
python ./launch.py --submit example.json start_time_scan
```

The fitting details could be defined in the `.json` configuration file.

After jobs are done, you should have the result files for each fitting in the `output_dir` defined in `.json` file.

Then you `hadd` and `delete` the files by
```
python ./launch.py --fetch example.json energy_bin
python ./launch.py --clean example.json energy_bin
```
You also can specify the output directory of *hadded* files by adding the path following `energy_bin`, the default is `./fetch`
