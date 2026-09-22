export PYTHONPATH=/insomnia001/depts/pmg/users/ic2465/MAG-E/
assemblers=$(python -m maggie query list --of assemblers)
targets=$(python -m maggie query list --of targets)

for t in $targets; do
  for a in $assemblers; do
    readarray -d '' aopts < <(python -m maggie query list --of options --within $a)
    for aopt in "${aopts[@]}"; do
        python -m maggie run --check-done --assembler $a --assembler-options $aopt --target $t
        if [ $? -eq 1 ]; then
          TARGET=$t AOPT=$aopt ASM=$a sbatch --export=ALL --job-name=$a assembly/run_assembler.sh
  	else
          echo done TARGET=$t AOPT=$aopt ASM=$a sbatch --export=ALL assembly/run_assembler.sh
	fi
    done
  done
done
