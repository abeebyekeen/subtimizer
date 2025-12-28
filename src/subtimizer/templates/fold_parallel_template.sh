#!/bin/bash
#SBATCH --job-name=AF-multi_fold_parallel
#SBATCH -p GPU4v100
#SBATCH --gres=gpu:2
#SBATCH -N 1
#SBATCH --mem 371552
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA (Adjust as needed)
# module load cuda118/toolkit/11.8.0 cuda118/blas/11.8.0 cuda118/fft/11.8.0
module load cuda/11.8.0

# Set ColabFold Path
# export PATH="/PathTo/colabfold/localcolabfold/colabfold-conda/bin:$PATH"

# Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_des

# Config
list_file="{file_path}"
start={start}
end={end}
max_parallel_jobs={max_parallel_jobs} # Matches GPU count usually

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

# Launch loop
line_num=0
while IFS= read -r des || [ -n "$des" ]; do
    line_num=$(( line_num + 1 ))
    
    # Filter lines
    if (( line_num < start )); then continue; fi
    if (( line_num > end )); then break; fi

    # Wait for slot
    while [ ${{#gpu_jobs[@]}} -ge $max_parallel_jobs ]; do
        check_completed_jobs
        sleep 10
    done

    # Launch
    gpu_id=$(find_available_gpu)
    
    echo "Launching $des on GPU $gpu_id..."
    
    # Sub-shell for the job
    (
        cd "$des"
        mkdir -p AFcomplex/top5complex
        round=5
        while [[ "$round" != 0 ]]; do
            mkdir -p AFcomplex/round_$round
            seed=$((1000 + RANDOM % 8999))
            fasta=$(ls *.fasta | head -n 1) # Assumes 1 output
            
            CUDA_VISIBLE_DEVICES=$gpu_id colabfold_batch \
                --num-recycle 10 --num-models 5 --random-seed "$seed" \
                --model-type auto --templates --amber --num-relax 3 --use-gpu-relax \
                "$fasta" AFcomplex/round_$round > /dev/null 2>&1
            
            # Copy result
            fasta_base=$(basename "$fasta" ".fasta")
            cp AFcomplex/round_$round/*_relaxed_rank_001_*.pdb AFcomplex/top5complex/ || true
            
            round=$(( round - 1 ))
        done
    ) &
    
    gpu_jobs[$gpu_id]=$!
    
    # Tiny sleep to ensure PID registered
    sleep 1

done < "$list_file"

wait
echo "All parallel folding jobs completed."
