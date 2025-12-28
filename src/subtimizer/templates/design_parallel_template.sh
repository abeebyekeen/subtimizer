#!/bin/bash
#SBATCH --job-name=mpnn_des_parallel
#SBATCH -p GPU4v100
#SBATCH --gres=gpu:4
#SBATCH -N 1
#SBATCH --mem 371552
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA
module load cuda/11.8.0

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mpnn_des

# Config matches input args
list_file="{file_path}"
config_file="{config_file}"
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

# Launch loop
line_num=0
while IFS= read -r des; do
    line_num=$(( line_num + 1 ))
    
    if (( line_num < start )); then continue; fi
    if (( line_num > end )); then break; fi

    # Wait for slot
    while [ ${{#gpu_jobs[@]}} -ge $max_parallel_jobs ]; do
        check_completed_jobs
        sleep 10
    done

    # Launch task
    gpu_id=$(find_available_gpu)
    echo "Launching $des on GPU $gpu_id..."
    
    # Lookup configuration from JSON using Python
    # This ensures "chains" and "fixed" are extracted correctly even if config file is complex
    # Fallback to defaults if JSON lookup fails or complex not found
    
    read chains_to_design fixed_positions <<< $(python -c "
import json, sys
try:
    with open('$config_file') as f:
        data = json.load(f)
        cfg = data.get('$des', {{}})
        print(cfg.get('chains_to_design', 'B'))
        print(cfg.get('fixed_positions', '4'))
except:
    print('B')
    print('4')
" | xargs)

    (
        cd "$des/AFcomplex/mpnn_des" || exit 1
        
        folder_with_pdbs="../top5complex/"
        output_dir="../mpnn_out"
        mkdir -p $output_dir

        path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
        path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
        path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
        
        python ${{MPNN_PATH}}/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains > /dev/null 2>&1
        python ${{MPNN_PATH}}/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design" > /dev/null 2>&1
        python ${{MPNN_PATH}}/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions" > /dev/null 2>&1

        CUDA_VISIBLE_DEVICES=$gpu_id python ${{MPNN_PATH}}/protein_mpnn_run.py \
            --jsonl_path $path_for_parsed_chains \
            --chain_id_jsonl $path_for_assigned_chains \
            --fixed_positions_jsonl $path_for_fixed_positions \
            --out_folder $output_dir \
            --num_seq_per_target 480 \
            --sampling_temp "0.1" \
            --batch_size 32 > "${{des}}_mpnn.log" 2>&1
    ) &

    gpu_jobs[$gpu_id]=$!
    sleep 1

done < "$list_file"

wait
echo "All parallel design jobs completed."
