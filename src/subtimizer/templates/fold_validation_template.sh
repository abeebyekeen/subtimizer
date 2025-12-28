#!/bin/bash
#SBATCH --job-name={complex_name}_val
#SBATCH -p GPUv100s
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --mem 360928
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load modules (Adjust as needed)
module load cuda/11.8.0

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_des

export CUDA_VISIBLE_DEVICES=0

echo "Running Validation Folding for {complex_name}"

# Navigate to the sequences directory
# Path relative to where script is submitted (the complex root)
WORK_DIR="AFcomplex/mpnn_out_clust_fold/seqs"
STRUCT_DIR="../structs"

if [ ! -d "$WORK_DIR" ]; then
    echo "Error: Directory $WORK_DIR not found."
    exit 1
fi

cd "$WORK_DIR"
mkdir -p "$STRUCT_DIR"

total_des=$(find . -maxdepth 1 -type d | wc -l)
current=0

# Loop through each design folder
for des_dir in */ ; do
    [ -d "$des_dir" ] || continue
    des_name=$(basename "$des_dir")
    
    pushd "$des_dir" > /dev/null
    
    # Skip if already done (check for result files)
    # Using heuristic from legacy script
    if ls *rank_001*.pdb 1> /dev/null 2>&1; then
       echo "Skipping $des_name (already done)"
       popd > /dev/null
       continue
    fi

    fasta=$(ls *.fasta | head -n 1)
    if [ -z "$fasta" ]; then
        echo "No fasta in $des_name"
        popd > /dev/null
        continue
    fi
    fasta_base=$(basename "$fasta" ".fasta")

    echo "Folding $des_name ($current / $total_des)..."

    # Run ColabFold
    # Less recycling/models than initial guess to save time (matches legacy settings usually)
    # Legacy 16: --num-recycle 2 --num-models 2
    colabfold_batch \
        --num-recycle 2 \
        --num-models 2 \
        --model-type auto \
        --templates \
        "$fasta" .
        
    # Copy top prediction to structs
    # Relative path from inside des_dir/ to structs/
    # (des_dir is in seqs/, structs/ is in mpnn_out_clust_fold/)
    # So path is ../../structs/
    cp ${{fasta_base}}_unrelaxed_rank_001_*.pdb ../../structs/ 2>/dev/null

    popd > /dev/null
    
    current=$((current + 1))
done

echo "Validation folding completed."
