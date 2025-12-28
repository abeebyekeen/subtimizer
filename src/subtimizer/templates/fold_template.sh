#!/bin/bash
#SBATCH --job-name={complex_name}_fold
#SBATCH -p GPUv100s
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --mem 360928
#SBATCH -t 3-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA if needed (Adjust version as per your HPC)
# module load cuda118/toolkit/11.8.0 cuda118/blas/11.8.0 cuda118/fft/11.8.0
module load cuda/11.8.0

# Set ColabFold Path (Adjust this path to your installation if it's not in PATH)
# export PATH="/PathTo/colabfold/localcolabfold/colabfold-conda/bin:$PATH"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_des

export CUDA_VISIBLE_DEVICES=0

echo "Launched colabfold for {complex_name}"
mkdir -p AFcomplex/top5complex

# Run 5 rounds of folding
for round in {{5..1}}; do
    echo "Starting Round $round..."
    mkdir -p AFcomplex/round_$round
    
    # Generate random seed
    seed=$((1000 + RANDOM % 8999))
    
    # Find fasta file
    fasta=$(ls *.fasta | head -n 1)
    if [ -z "$fasta" ]; then
        echo "Error: No fasta file found in $(pwd)"
        exit 1
    fi
    fasta_base=$(basename "${{fasta}}" ".fasta")

    colabfold_batch \
        --num-recycle 10 \
        --num-models 5 \
        --random-seed "$seed" \
        --model-type auto \
        --templates \
        --amber \
        --num-relax 3 \
        --use-gpu-relax \
        "${{fasta}}" AFcomplex/round_$round

    echo "Round $(( 6 - round )) completed."
    
    # Copy top prediction
    cp AFcomplex/round_$round/${{fasta_base}}_relaxed_rank_001_*.pdb AFcomplex/top5complex/ || echo "Warning: Could not copy ranking file."
done
