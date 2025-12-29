import os
import time
import pkgutil
import csv
from subtimizer.manager import JobManager

def run_design(file_path: str, max_jobs: int = 4, start: int = 1, end: int = None, mode: str = 'batch'):
    """
    Run ProteinMPNN design process
    """
    
    with open(file_path, 'r') as f:
        all_complexes = [line.strip() for line in f if line.strip()]
        
    if end is None:
        end = len(all_complexes)
        
    complexes = all_complexes[start-1 : end]

    # Load configuration
    import json
    config_file = "design_config.json"
    config_abs_path = os.path.abspath(config_file)
    configs = {}
    if os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                configs = json.load(f)
        except Exception as e:
            print(f"Warning: Error reading {config_file}: {e}")

    if mode == 'parallel':
        print(f"Submitting PARALLEL design job for {len(complexes)} complexes (Index {start} to {end})...")
        _submit_parallel_job(file_path, max_jobs, start, end, config_abs_path)
        return

    manager = JobManager(max_jobs=max_jobs)
    
    for complex_name in complexes:
        manager.wait_for_slot()
        
        mpnn_folder = os.path.join(complex_name, "AFcomplex", "mpnn_des")
        if not os.path.exists(mpnn_folder):
            print(f"Error: {mpnn_folder} does not exist. Run 'setup --type mpnn' first.")
            continue
            
        script_path = os.path.join(mpnn_folder, "run_mpnn.sh")
        
        # Get config for this complex (or default)
        cfg = configs.get(complex_name, {'chains_to_design': 'B', 'fixed_positions': '4'})
        
        if _write_mpnn_script(script_path, complex_name, cfg['chains_to_design'], cfg['fixed_positions']):
            wd = os.getcwd()
            os.chdir(mpnn_folder)
            manager.submit_job("run_mpnn.sh", job_name=f"mpnn_{complex_name}")
            os.chdir(wd)
            time.sleep(2)
        else:
            print(f"Skipping {complex_name} due to script generation error.")

def _write_mpnn_script(path, complex_name, chains, fixed):
    try:
        from subtimizer.utils import get_template_content
        template_content = get_template_content('design_template.sh')
    except Exception as e:
        print(f"Error loading template: {e}")
        return False

    content = template_content.format(
        complex_name=complex_name, 
        chains_to_design=chains, 
        fixed_positions=fixed
    )
    with open(path, 'w') as f:
        f.write(content)
    return True

def _submit_parallel_job(file_path, max_parallel_jobs, start, end, config_file_path):
    from subtimizer.utils import get_template_content
    try:
        template_content = get_template_content('design_parallel_template.sh')
    except Exception as e:
        print(f"Error loading parallel template: {e}")
        return

    script_name = "run_design_parallel.sh"
    content = template_content.format(
        file_path=os.path.abspath(file_path),
        start=start,
        end=end,
        max_parallel_jobs=max_parallel_jobs,
        config_file=config_file_path
    )
    
    with open(script_name, 'w') as f:
        f.write(content)
    
    print(f"Created {script_name}. Submitting...")
    job_id = JobManager.submit_slurm_job(script_name)
    print(f"Submitted parallel design job ID: {job_id}")
