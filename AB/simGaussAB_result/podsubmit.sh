#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=4
#SBATCH --time=700:00:00
#SBATCH --job-name=gaussAB
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

cd $SLURM_SUBMIT_DIR
/bin/hostname
OPENMM_CPU_THREADS=4
#export PATH=/home/mnguyen/bin/lammps/lammps-12Dec18/bin/:$PATH
export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
python main.py

