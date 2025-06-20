#!/bin/bash
#SBATCH -J ms_causal
#SBATCH -o ms_expcausal_output.txt
#SBATCH -e ms_causal_errors.txt
#SBATCH --partition=sapphire
#SBATCH -t 24:55:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=keithbarnatchez@g.harvard.edu
#SBATCH --mail-type=ALL

module load R

Rscript sim-main-dr.R
