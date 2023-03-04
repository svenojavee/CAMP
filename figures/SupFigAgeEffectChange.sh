#!/bin/bash

#SBATCH --job-name=sFig
#SBATCH --partition=defaultp
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 150gb
#SBATCH --time 0-10:00:00
#SBATCH --output=/nfs/scistore13/robingrp/sojavee/figureData/SupFigAgeEffect.log

module add R

Rscript /nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/SupFigAgeEffectChange.R
