# Omega-a fitter for g-2 analysis

## Introduction
Perform omega-a fitting on the wiggle histogram. This fitter have features:

1. Integrated Omega-a `Blinder`
2. Easily submit fitting jobs to `inpac condor` for seed scan, start time scan, adhoc scan, etc.


To use this fitter, prepare the ROOT files with

1. Wiggle histogram 
1. Lost muon integral histogram
1. Initial values saved as `TArrayD` array or `TF1` function

Then write up the configuration `json` file with necessary fitting setup and submit the fitter to `condor (inpac)` or run `local`


## An example with inpac-condor

Clone git and build with Run4 `blinding`

```
mkdir Fitter_wa && cd Fitter_wa
git clone https://github.com/heymanwasup/OmegaFitting.git
source ./OmegaFitting/initialize_inpac.sh run4
```

The `initialize_inpac.sh` script will build the project as well as make run directory and link some configuration file into it. If you already have those files and just want to rebuild the project, you can run
```
source ./OmegaFitting/initialize_inpac.sh run4 false
```
in the Fitter_wa directory.

Now, you can find an example of configuration file as `example.json`, you can perform several tasks with this configuration file:

1. Submit seed scan:
```
python ./launch.py --submit example.json run2_seed_scan
```
This command will luanch `200 seeds * T-/A-method * Run2all/Run2C = 800` jobs in total with `inpac condor` with 28 parameters model

2. Submit calorimeter scan
```
python ./launch.py --submit example.json run2_calorimeter_scan
```

The fitting details are defined in the `example.json` file.

3. (other tasks to be added to the example.json)
