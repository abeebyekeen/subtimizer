#!/bin/bash
#SBATCH --job-name={complex_name}_mpnn
#SBATCH -p GPUv100s
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --mem 360928
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA (Adjust as needed)
# module load cuda112/toolkit/11.2.0 cuda112/blas/11.2.0 cuda112/fft/11.2.0
module load cuda/11.8.0

export CUDA_VISIBLE_DEVICES=0

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mpnn_des

# Set paths
folder_with_pdbs="../top5complex/"
output_dir="../mpnn_out"
mkdir -p $output_dir

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"

# Configuration injected by subtimizer
chains_to_design="{chains_to_design}"
fixed_positions="{fixed_positions}"

echo "Configuration: Chains=$chains_to_design, Fixed=$fixed_positions"

# Ensure MPNN_PATH is set
if [ -z "$MPNN_PATH" ]; then
    echo "Error: MPNN_PATH is not set. Please set it in your environment or this script."
    exit 1
fi

echo "Running helper scripts..."
python ${{MPNN_PATH}}/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python ${{MPNN_PATH}}/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python ${{MPNN_PATH}}/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

echo "Running ProteinMPNN..."
python ${{MPNN_PATH}}/protein_mpnn_run.py \
    --jsonl_path $path_for_parsed_chains \
    --chain_id_jsonl $path_for_assigned_chains \
    --fixed_positions_jsonl $path_for_fixed_positions \
    --out_folder $output_dir \
    --num_seq_per_target 480 \
    --sampling_temp "0.1" \
    --batch_size 32

echo "Completed design for {complex_name}"
