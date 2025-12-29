import os
import shutil
import time
from subtimizer.manager import JobManager

def run_clustering(file_path: str, max_jobs: int = 4):
    """
    Run CD-HIT clustering
    """
    print(f"Running filtering/clustering for complexes in {file_path}")
    
    with open(file_path, 'r') as f:
        complexes = [line.strip() for line in f if line.strip()]

    manager = JobManager(max_jobs=max_jobs)
    
    for complex_name in complexes:
        manager.wait_for_slot()
        
        # Paths
        af_complex_dir = os.path.join(complex_name, "AFcomplex")
        mpnn_out_dir = os.path.join(af_complex_dir, "mpnn_out", "seqs")
        clust_dir = os.path.join(af_complex_dir, "mpnn_out_clust")
        
        source_fasta = os.path.join(mpnn_out_dir, "all_design.fa")
        
        if not os.path.exists(source_fasta):
            print(f"Skipping {complex_name}: Input fasta {source_fasta} not found.")
            continue
            
        # Setup directories
        if os.path.exists(clust_dir):
            shutil.rmtree(clust_dir)
        os.makedirs(clust_dir, exist_ok=True)
        
        shutil.copy(source_fasta, os.path.join(clust_dir, "all_design.fa"))
        
        # Prepare script
        script_path = os.path.join(clust_dir, "run_cdhit.sh")
        _write_cdhit_script(script_path, complex_name)
        
        # Submit
        wd = os.getcwd()
        os.chdir(clust_dir)
        manager.submit_job("run_cdhit.sh", job_name=f"cdhit_{complex_name}")
        os.chdir(wd)
        
        time.sleep(1)

    print("Waiting for clustering jobs to finish before summarizing...")
    pass

def summarize_clusters(file_path: str):
    """
    Waits for clustering jobs to finish and generates summary.
    """
    with open(file_path, 'r') as f:
        complexes = [line.strip() for line in f if line.strip()]
    
    # timeout: 5 minutes per complex, maxed at 30 minutes
    num_complexes = len(complexes)
    calculated_timeout = num_complexes * 5 * 60
    max_timeout = 30 * 60
    timeout = min(calculated_timeout, max_timeout)
    
    print(f"Waiting for clustering jobs to complete (Timeout: {timeout/60:.1f}m)...")
    
    start_time = time.time()
    
    while True:
        pending = []
        for folder in complexes:
            # Check for output file (success) and log presence
            clust_file = os.path.join(folder, "AFcomplex", "mpnn_out_clust", "all_design_clustered.fa")
            log_file = os.path.join(folder, "AFcomplex", "mpnn_out_clust", "cdhit.log")
            
            is_done = False
            if os.path.exists(log_file):
                # check for completion in log
                try:
                    with open(log_file, 'r') as lf:
                        if "finished" in lf.read():
                            is_done = True
                except:
                    pass
            
            if not is_done:
                pending.append(folder)
        
        if not pending:
            print("All clustering jobs completed.")
            break
            
        elapsed = time.time() - start_time
        if elapsed > timeout:
            print("\nTimeout reached! Proceeding with summary for completed jobs only.")
            print(f"Pending/Incomplete: {len(pending)} complexes ({', '.join(pending[:5])}...)")
            break
            
        # Update user every check
        print(f"Waiting for {len(pending)}/{len(complexes)} jobs... ({int(elapsed)}s elapsed)", end='\r')
        time.sleep(5)

    print("\nSummarizing clustering results...")
    
    summary_file = "cluster_summary.dat"
    with open(summary_file, "w") as out:
        for folder in complexes:
            log_file = os.path.join(folder, "AFcomplex", "mpnn_out_clust", "cdhit.log")
            if os.path.exists(log_file):
                try:
                    with open(log_file) as lf:
                        for line in lf:
                            if "finished" in line:
                                parts = line.strip().split()
                                if len(parts) >= 4:
                                    out.write(f"{folder} :   {parts[2]}  {parts[3]}\n")
                except Exception as e:
                    print(f"Error reading log for {folder}: {e}")
            else:
                 pass
    
    print(f"Summary written to {summary_file}")

import pkgutil

def _write_cdhit_script(path, complex_name):
    try:
        template_bytes = pkgutil.get_data('subtimizer.templates', 'cluster_template.sh')
        if template_bytes is None:
             raise FileNotFoundError("Template cluster_template.sh not found.")
        template_content = template_bytes.decode('utf-8')
    except Exception as e:
        print(f"Error loading template: {e}")
        return False

    content = template_content.format(complex_name=complex_name)
    with open(path, 'w') as f:
        f.write(content)
    return True
