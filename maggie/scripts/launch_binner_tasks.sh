# Map each assembly using each of the mapping tools 
# run with all supplied options over all samples.

# List of all target samples
targets=$(python -m maggie query list --of targets)

# List of assembler tools
assemblers=$(python -m maggie query list --of assemblers)

# List of all mapper tools
mappers=$(python -m maggie query list --of mappers)

# List of binning modes. For modes other than single, 
# we may need to map other samples to the target assemblies
modes=$(python -m maggie query list --of modes)

# List binners
binners=$(python -m maggie query list --of binners)

# Loop through each mapping task (product over assemblers, mappers, and their options).
for a in MEGAHIT; do 
  readarray -d '' aopts < <(python -m maggie query list --of options --within $a)
  for aopt in "${aopts[@]}"; do 
    for m in $mappers; do 
      readarray -d '' mopts < <(python -m maggie query list --of options --within $m)
      for mopt in "${mopts[@]}"; do 
        for b in VAMB; do
          readarray -d '' bopts < <(python -m maggie query list --of options --within $b)
          for bopt in "${bopts[@]}"; do
            for t in $targets; do 
              for d in $modes; do
                python -m maggie run --check-done --assembler $a --assembler-options $aopt --target $t --mapper $m --mapper-options $mopt --binner $b --binner-options $bopt --binning-mode $d
                if [ $? -eq 1 ]; then
		              if [ $b == "COMEBin" ]; then
                    echo ASM=$a AOPT=$aopt MAP=$m MOPT=$mopt BOPT=$bopt TARGET=$t MODE=$d BIN=$b sbatch --export=ALL binning/run_binner_gpu.sh
                    ASM=$a AOPT=$aopt MAP=$m MOPT=$mopt BOPT=$bopt TARGET=$t MODE=$d BIN=$b sbatch --export=ALL,MAMBA_NO_LOCK=1,CONDA_NO_LOCK=1 binning/run_binner_gpu.sh
		              else
                    echo ASM=$a AOPT=$aopt MAP=$m MOPT=$mopt BOPT=$bopt TARGET=$t MODE=$d BIN=$b sbatch --export=ALL binning/run_binner.sh
                    ASM=$a AOPT=$aopt MAP=$m MOPT=$mopt BOPT=$bopt TARGET=$t MODE=$d BIN=$b sbatch --export=ALL binning/run_binner.sh
	                fi
                else
                  echo $a $aopt $m $mopt $b $bopt $m $t done
                fi
              done
            done
          done
        done
      done
    done
  done
done
