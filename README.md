# Figure Reproduction Scripts

This repository contains two types of scripts to reproduce figures. They have been tested on R version 4.4.2.
on a x86_64-pc-linux-gnu machine.

## 1. run_*-scripts
Direct execution scripts that include all calculations. The user should be aware that this might take time
depending on the simulation setting.

Default simulation settings are estimated with n_simulation = 100 and can be changed in the scripts.

```r
source("scripts/run_figure_*.R")
```


## 2. targets_*-scripts
Scripts using the targets framework. The scripts use the target package to load prebuild binary objects 
to reproduce the figures of the master thesis. These scripts are guaranteed be very fast compared to the 
run_*-scripts.

Download and extract the targets folder into the project directory from the following link: 

```r
source("scripts/targets_figure_*.R")
```

## Output
Both script types generate figures in scripts/figures/

## Dependencies to run the scripts
ggplot2
gridExtra
grid
compositions
robCompositions
zCompositions
targets
dplyr

Running the scripts will automatically install missing packages if they are not installed on the users R version.