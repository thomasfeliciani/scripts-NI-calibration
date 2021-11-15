#!/bin/bash

#SBATCH --job-name=NIcalibration
#SBATCH --time=7-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3000
#SBATCH --mem=3GB

ml R/3.4.4-foss-2018a-X11-20180131
Rscript metaHandler_1c.R
