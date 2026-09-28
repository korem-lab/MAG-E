#!/bin/bash 
#SBATCH --job-name=mapper
#SBATCH --time=6:0:0
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --nodelist=ins0[80-92]
#SBATCH --cpus-per-task=8
#SBATCH --account pmg

if [ $STAGE == 'prep' ]; then
    python -m maggie run --assembler $ASM --assembler-options $AOPT --target $TARGET --threads 8 --coverage $MAP --coverage-options $MOPT --stage $STAGE
else
    python -m maggie run --assembler $ASM --assembler-options $AOPT --target $TARGET --threads 8 --coverage $MAP --coverage-options $MOPT --cov-sample $SAMPLE --stage $STAGE
fi
