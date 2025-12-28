#!/bin/bash
#SBATCH --job-name=AF-val_fold_parallel
#SBATCH -p GPU4v100
#SBATCH --gres=gpu:4
#SBATCH -N 1
#SBATCH --mem 371552
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA (Adjust as needed)
module load cuda/11.8.0

# Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_des

# Config
list_file="{file_path}"
start={start}
end={end}
max_parallel_jobs={max_parallel_jobs}

declare -A gpu_jobs=()

# Find available GPU
find_available_gpu() {{
    for gpu in $(seq 0 $((max_parallel_jobs - 1))); do
        if [[ -z ${{gpu_jobs[$gpu]}} ]]; then
            echo $gpu
            return
        fi
    done
}}

# Check for completed jobs
check_completed_jobs() {{
    for gpu in "${{!gpu_jobs[@]}}"; do
        local pid=${{gpu_jobs[$gpu]}}
        if ! kill -0 $pid 2>/dev/null; then
            unset gpu_jobs[$gpu]
        fi
    done
}}

line_num=0
while IFS= read -r complex || [ -n "$complex" ]; do
    line_num=$(( line_num + 1 ))
    
    # Filter lines
    if (( line_num < start )); then continue; fi
    if (( line_num > end )); then break; fi
    
    echo "Processing complex: $complex"
    
    # Validation folding needs to run on all designs within the complex
    WORK_DIR="$complex/AFcomplex/mpnn_out_clust_fold/seqs"
    STRUCT_DIR="../structs"
    
    if [ ! -d "$WORK_DIR" ]; then
        echo "Warning: $WORK_DIR not found. Skipping."
        continue
    fi
    
    # Find all design directories
    # Flatten loop: process designs one by one using the GPU pool
    for des_dir in "$WORK_DIR"/*/ ; do
        [ -d "$des_dir" ] || continue
        
        # Wait for slot
        while [ ${{#gpu_jobs[@]}} -ge $max_parallel_jobs ]; do
            check_completed_jobs
            sleep 5
        done
        
        gpu_id=$(find_available_gpu)
        
        # Launch Job
        (
            des_name=$(basename "$des_dir")
            # Navigate to design dir
            cd "$des_dir"
            
            # Skip if done (heuristic)
            if ls *rank_001*.pdb 1> /dev/null 2>&1; then
                 # echo "Skipping $complex/$des_name (done)"
                 exit 0
            fi

            fasta=$(ls *.fasta | head -n 1)
            if [ -z "$fasta" ]; then
                 exit 0
            fi
            fasta_base=$(basename "$fasta" ".fasta")
            
            # Ensure struct dir exists (relative to design dir: ../../structs)
            mkdir -p ../../structs/

            # Run ColabFold (Validation Settings)
            CUDA_VISIBLE_DEVICES=$gpu_id colabfold_batch \
                --num-recycle 2 --num-models 2 \
                --model-type auto --templates \
                "$fasta" . > /dev/null 2>&1
            
            # Copy result
            # Relative path from des_dir (inside seqs/) to ../../structs (inside mpnn_out_clust_fold/)
            cp ${{fasta_base}}_unrelaxed_rank_001_*.pdb ../../structs/ 2>/dev/null
            
            echo "Finished $complex/$des_name on GPU $gpu_id"
        ) &
        
        gpu_jobs[$gpu_id]=$!
        sleep 0.5
    done


done < "$list_file"

wait
echo "All parallel AlphaFold jobs completed."
