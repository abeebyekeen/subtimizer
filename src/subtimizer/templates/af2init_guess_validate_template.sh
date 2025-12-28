#!/bin/bash
#SBATCH --job-name={complex_name}_valid
#SBATCH --nodes=1
#SBATCH -p GPU4v100
#SBATCH --gres=gpu:4
#SBATCH --mem 371552
#SBATCH -t 6-23:45:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Load CUDA (Adjust as needed)
module load cuda/11.8.0
export CUDA_VISIBLE_DEVICES=0

source $(conda info --base)/etc/profile.d/conda.sh

conda activate af2_des


echo "Running Validation for {complex_name}"

# NOTE: The implementation assumes user has configured dl_binder_design
PREDICT_PY="{dl_binder_path}"

if [ ! -f "$PREDICT_PY" ]; then
    echo "Error: predict.py not found at $PREDICT_PY"
    echo "Please edit the script or configuration to point to dl_binder_design."
    exit 1
fi

python "$PREDICT_PY" \
    -pdbdir ../af2_init_guess_in \
    -outpdbdir ../af2_init_guess_out.rec8 \
    -scorefilename af2score.dat \
    -recycle 8
