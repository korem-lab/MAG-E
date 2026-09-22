# List of all target samples
targets=$(python -m maggie query list --of targets)

# List of assembler tools
assemblers=$(python -m maggie query list --of assemblers)

# List of all mapper tools
mappers=$(python -m maggie query list --of mappers)

for a in $assemblers; do 
  readarray -d '' aopts < <(python -m maggie query list --of options --within $a)
  for aopt in "${aopts[@]}"; do 
    for m in $mappers; do 
      readarray -d '' mopts < <(python -m maggie query list --of options --within $m)
      for mopt in "${mopts[@]}"; do 
        for t in $targets; do 
          python -m maggie run --check-done --assembler $a --assembler-options $aopt --target $t --mapper $m --mapper-options $mopt
          if [ $? -eq 1 ]; then 
            TARGET=$t ASM=$a AOPT=$aopt MOPT=$mopt MAP=$m STAGE=prep sbatch --export=ALL mapping/run_mapper.sh
          fi 
            exit 0
        done
      done
    done
  done
done