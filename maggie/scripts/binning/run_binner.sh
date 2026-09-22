#!/bin/bash 
#SBATCH --job-name=binner
#SBATCH --time=12:0:0
#SBATCH --mem=64G
#SBATCH --exclude=ins006,ins005
#SBATCH --cpus-per-task=8
#SBATCH --account pmg

python -m maggie run --assembler $ASM --assembler-options $AOPT --mapper $MAP --mapper-options $MOPT --binner $BIN --binner-options $BOPT --binning-mode $MODE --target $TARGET --threads 8
