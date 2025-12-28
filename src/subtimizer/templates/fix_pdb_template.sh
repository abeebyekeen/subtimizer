#!/bin/bash
#SBATCH --job-name={complex_name}_fix
#SBATCH -p 512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem 501760  
#SBATCH -t 04:00:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

echo "Fixing PDBs for {complex_name}"

subtimizer internal-fix-pdb --dir {target_dir}
