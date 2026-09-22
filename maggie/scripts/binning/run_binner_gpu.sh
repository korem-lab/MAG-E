#!/bin/bash 
#SBATCH --job-name=binner
#SBATCH --time=4:0:0
#SBATCH --mem=128G
#SBATCH --nodelist=ins094
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --account pmg

python -m maggie run --assembler $ASM --assembler-options $AOPT --mapper $MAP --mapper-options $MOPT --binner $BIN --binner-options $BOPT --binning-mode $MODE --target $TARGET --threads 16 --stage main
