import os
import shutil
import glob
from subtimizer.manager import JobManager

import pkgutil

def run_validation(file_path: str, max_jobs: int = 4, binder_path: str = None, start: int = 1, end: int = None):
    """
    Run AF2 Initial Guess validation
    """
    print(f"Running AF2 Init Guess for complexes in {file_path}")
    
    # Resolve binder path
    if binder_path is None:
        binder_path = "path_to/dl_binder_design/af2_initial_guess/predict.py"
        print("Warning: --binder-path not provided. Using default placeholder.")
    
    with open(file_path, 'r') as f:
        all_complexes = [line.strip() for line in f if line.strip()]
        
    # Filter by start/end (1-based)
    if end is None:
        end = len(all_complexes)
    
    complexes = all_complexes[start-1 : end]
    print(f"Processing {len(complexes)} complexes (Index {start} to {end})...")

    manager = JobManager(max_jobs=max_jobs)
    
    for complex_name in complexes:
        manager.wait_for_slot()
        
        # Determine Mode
        standard_fold_dir = os.path.join(complex_name, "AFcomplex", "mpnn_out_clust_fold")
        flat_top5_dir = os.path.join(complex_name, "top5complex")
        
        mode = "unknown"
        if os.path.exists(standard_fold_dir):
            mode = "standard"
            work_dir = standard_fold_dir
        elif os.path.exists(flat_top5_dir):
            mode = "flat"
            work_dir = complex_name
        else:
            print(f"Skipping {complex_name}: No valid directory found (Standard or Flat).")
            continue

        # Paths relative to work_dir
        init_guess_dir = os.path.join(work_dir, "af2_init_guess.rec8")
        init_guess_out = os.path.join(work_dir, "af2_init_guess_out.rec8")
        
        # Setup directories
        if os.path.exists(init_guess_dir):
            shutil.rmtree(init_guess_dir)
        os.makedirs(init_guess_dir, exist_ok=True)
        
        if os.path.exists(init_guess_out):
            shutil.rmtree(init_guess_out)
        os.makedirs(init_guess_out, exist_ok=True)

        # Prepare worker script
        script_path = os.path.join(init_guess_dir, "run_validate.sh")
        # Template relies on relative paths (../af2_init_guess_in), so structure must hold
        if _write_validation_script(script_path, complex_name, work_dir, binder_path):
            wd = os.getcwd()
            os.chdir(init_guess_dir)
            manager.submit_job("run_validate.sh", job_name=f"valid_{complex_name}")
            os.chdir(wd)

def _write_validation_script(path, complex_name, fold_dir, dl_binder_path):
    """
    Writes validation script using external template
    """
    try:
        template_bytes = pkgutil.get_data('subtimizer.templates', 'af2init_guess_validate_template.sh')
        if template_bytes is None:
             raise FileNotFoundError("af2init_guess_validate_template.sh not found")
        template_content = template_bytes.decode('utf-8')
    except Exception as e:
        print(f"Error loading template: {e}")
        return False

    content = template_content.format(
        complex_name=complex_name,
        dl_binder_path=dl_binder_path
    )
    
    with open(path, 'w') as f:
        f.write(content)
    return True
