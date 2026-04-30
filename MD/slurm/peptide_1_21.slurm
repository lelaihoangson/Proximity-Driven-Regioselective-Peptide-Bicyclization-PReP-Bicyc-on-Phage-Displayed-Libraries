#!/bin/bash
#SBATCH --job-name=MD-run
#SBATCH --output=peptide_1_21.slurm.oe%j
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # 1 MPI rank
#SBATCH --cpus-per-task=12           # 12 OpenMP threads
#SBATCH --gres=gpu:1                 # 1 GPU
#SBATCH --partition=gpu              # Use the GPU partition
#SBATCH --mem=4096M                  # Memory per node

##SBATCH --mail-type=ALL
##SBATCH --mail-user=lelaihoangson@tamu.edu

# Load modules
ml purge
ml GCC/13.2.0 
ml OpenMPI/4.1.6
ml GROMACS/2024.4-CUDA-12.4.1

export OMP_NUM_THREADS=12

# Move to working directory
cd /scratch/user/lelaihoangson/PD1/MD/peptide_1_21

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10 -n index.ndx
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 12

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 10 -n index.ndx
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 12

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_100.tpr -maxwarn 10 -n index.ndx
gmx mdrun -deffnm md_0_100 -ntmpi 1 -ntomp 12
