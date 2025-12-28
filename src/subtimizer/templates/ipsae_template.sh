#!/bin/bash
#SBATCH --job-name=ipsae_run
#SBATCH -p 512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task={max_jobs}
#SBATCH --mem 501760
#SBATCH -t 04:00:00
#SBATCH -o ipsae_job_%j.out
#SBATCH -e ipsae_job_%j.err

echo "Running ipSAE evaluation..."
echo "List: {list_file}"
echo "Jobs: {max_jobs}"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate subtimizer_env

subtimizer internal-ipsae --file "{list_file}" --pae-cutoff {pae_cutoff} --dist-cutoff {dist_cutoff} --max-jobs {max_jobs} --start {start} --end {end}
