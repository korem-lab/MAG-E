#!/bin/bash 
#SBATCH --job-name=mapper
#SBATCH --time=6:0:0
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg

if [ $STAGE == 'prep' ]; then
    python -m maggie run --assembler $ASM --assembler-options $AOPT --target $TARGET --threads 8 --mapper $MAP --mapper-options $MOPT --stage $STAGE
else
    python -m maggie run --assembler $ASM --assembler-options $AOPT --target $TARGET --threads 8 --mapper $MAP --mapper-options $MOPT --map-sample $SAMPLE --stage $STAGE
fi
