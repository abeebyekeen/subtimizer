import os
import time
import pkgutil
from subtimizer.manager import JobManager

def run_folding(file_path: str, max_jobs: int = 4, start: int = 1, end: int = None, mode: str = 'batch', stage: str = 'initial'):
    """
    Execute AF-multimer/AF2-initial-guess
    mode: 'batch' (one job per complex) or 'parallel' (one job for all, multi-gpu).
    stage: 'initial' (generating initial structure) or 'validation' (AF2-initial-guess)
    """
    
    with open(file_path, 'r') as f:
        all_complexes = [line.strip() for line in f if line.strip()]

    # Filter by start/end (1-based indexing)
    if end is None:
        end = len(all_complexes)
    
    # Adjust for 0-based indexing
    complexes = all_complexes[start-1 : end]
    
    if mode == 'parallel':
        # Parallel Mode: Submit ONE job that handles the pool internally
        print(f"Submitting PARALLEL job for {len(complexes)} complexes (Index {start} to {end})...")
        if stage == 'validation':
             _submit_parallel_job(file_path, max_jobs, start, end, template_name="fold_validation_parallel_template.sh")
        else:
             _submit_parallel_job(file_path, max_jobs, start, end, template_name="fold_parallel_template.sh")
        return

    # Batch Mode: One job per complex
    manager = JobManager(max_jobs=max_jobs)
    
    for complex_name in complexes:
        # Check if already done
        if stage == 'initial':
            done_file = os.path.join(complex_name, "AFcomplex", "round_1", f"{complex_name}.done.txt") # heuristic
        else:
            done_file = os.path.join(complex_name, "AFcomplex", "mpnn_out_clust_fold", "validation.done")

        if os.path.exists(done_file):
            print(f"{complex_name} already completed. Skipping.")
            continue
            
        manager.wait_for_slot()
        
        # Determine script name and template
        if stage == 'validation':
            worker_script_path = os.path.join(complex_name, "run_validation_fold.sh")
            template_name = "fold_validation_template.sh"
            job_name = f"val_{complex_name}"
        else:
            worker_script_path = os.path.join(complex_name, "run_fold.sh")
            template_name = "fold_template.sh"
            job_name = f"fold_{complex_name}"

        if _write_worker_script(worker_script_path, complex_name, template_name):
            # Submit only if script writing succeeded
            wd = os.getcwd()
            os.chdir(complex_name)
            script_name = os.path.basename(worker_script_path)
            print(f"Submitting: {script_name} ({stage})")
            manager.submit_job(script_name, job_name=job_name)
            os.chdir(wd)
            
            # Small delay to ensure SLURM registers the job
            time.sleep(2)
        else:
            print(f"Skipping {complex_name} due to script generation error.")

def _write_worker_script(path, complex_name, template_name="fold_template.sh"):
    """
    Writes the specific folding script for a complex using the package template.
    """
    try:
        from subtimizer.utils import get_template_content
        template_content = get_template_content(template_name)
    except Exception as e:
        print(f"Error loading template: {e}")
        return False

    # Replace placeholders
    content = template_content.format(complex_name=complex_name)
    
    with open(path, 'w') as f:
        f.write(content)
    return True

def _submit_parallel_job(file_path, max_parallel_jobs, start, end, template_name="fold_parallel_template.sh"):
    """
    Creates and submits a single script that runs multiple GPU tasks in parallel.
    Uses fold_parallel_template.sh or fold_validation_parallel_template.sh
    """
    try:
    # Robust template loading
    from subtimizer.utils import get_template_content
    try:
        template_content = get_template_content(template_name)
    except Exception as e:
        print(f"Error loading parallel template: {e}")
        return

    # Create the single master script in the current directory
    # Distinguish script names and make unique by range
    if "validation" in template_name:
        script_name = f"run_validation_parallel_{start}_{end}.sh"
    else:
        script_name = f"run_fold_parallel_{start}_{end}.sh"
        
    content = template_content.format(
        file_path=os.path.abspath(file_path),
        start=start,
        end=end,
        max_parallel_jobs=max_parallel_jobs
    )
    
    with open(script_name, 'w') as f:
        f.write(content)
    
    print(f"Created {script_name}. Submitting...")
    job_id = JobManager.submit_slurm_job(script_name)
    print(f"Submitted parallel job ID: {job_id}")
