script_qe = """#!/bin/sh
#$ -cwd
#$ -V -S /bin/bash
####$ -V -S /bin/zsh
#$ -N Si8_CrySPY_ID
#$ -pe smp 20
####$ -q ibis1.q
####$ -q ibis2.q

mpirun -np {mpi} pw.x < {prefix}.in > {prefix}.out


if [ -e "CRASH" ]; then
    sed -i -e '3 s/^.*$/skip/' stat_job
    exit 1
fi

sed -i -e '3 s/^.*$/done/' stat_job
"""

script_mlp = """#!/bin/sh

# ---------- ASE
python ase_in.py

# ---------- CrySPY
sed -i -e '3 s/^.*$/done/' stat_job
"""
