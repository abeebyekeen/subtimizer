#!/bin/bash
#SBATCH --job-name={complex_name}_cdhit
#SBATCH -p 512GB
#SBATCH -N 1
#SBATCH --mem 501760
#SBATCH -t 23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CD-HIT module (Adjust as needed for your HPC)
module load cdhit/4.8.1

echo "Running CD-HIT for {complex_name}"

# CD-HIT command
# -i: Input file
# -o: Output file
# -c: Sequence identity threshold (1.0 = 100%)
# -M: Memory limit (in MB)
# -T: Number of threads (0 = all available)
# -l: Length of throw_away_sequences

cd-hit -i all_design.fa -o all_design_clustered.fa -c 1.0 -M 32000 -T 0 -l 5 > cdhit.log
