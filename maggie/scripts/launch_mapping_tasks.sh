#!/bin/bash 
#SBATCH --job-name=launch
#SBATCH --time=6:0:0
#SBATCH --mem=4G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=2
export PYTHONPATH=/insomnia001/depts/pmg/users/ic2465/MAG-E


# Map each assembly using each of the mapping tools 
# run with all supplied options over all samples.

# List of all target samples
targets=$(python -m maggie query list --of targets)

# List of assembler tools
assemblers=$(python -m maggie query list --of assemblers)

# List of all mapper tools
mappers=$(python -m maggie query list --of mappers)

# Loop through each mapping task (product over assemblers, mappers, and their options).
for a in metaSPAdes; do 
  readarray -d '' aopts < <(python -m maggie query list --of options --within $a)
  for aopt in "${aopts[@]}"; do 
    for m in $mappers; do 
      readarray -d '' mopts < <(python -m maggie query list --of options --within $m)
        for mopt in "${mopts[@]}"; do 
          for t in $targets; do 
            if [[ $t == SRR12344455 ]]; then continue; fi
            for mode in all; do
              # For each mode get the list of samples used to bin the target.
              readarray -d '' samples < <(python -m maggie query list --of samples --within $mode --target $t)
              for s in $samples; do 
                #python -m maggie run --check-done --stage main --assembler $a --assembler-options $aopt --target $t --map-sample $s --mapper $m --mapper-options $mopt
                #if [ $? -eq 1 ]; then 
                  TARGET=$t SAMPLE=$s ASM=$a AOPT=$aopt MOPT=$mopt MAP=$m STAGE=main sbatch --export=ALL --job-name $m mapping/run_mapper.sh
		  sleep 0.1
                #fi 
            done
          done
        done
      done
    done
  done
done
