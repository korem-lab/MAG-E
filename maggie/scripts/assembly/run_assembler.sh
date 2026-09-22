#!/bin/bash 
#SBATCH --job-name=mh
#SBATCH --time=24:0:0
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --account pmg

python -m maggie run --assembler $ASM --assembler-options $AOPT --target $TARGET --threads 32
