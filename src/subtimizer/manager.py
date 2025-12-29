import subprocess
import time
import os

class JobManager:
    """
    Manages SLURM job submissions and monitoring.
    """
    def __init__(self, max_jobs=4, user=None):
        self.max_jobs = max_jobs
        self.user = user or os.getenv('USER')
        
    def get_running_jobs(self):
        """Returns the number of running jobs for the user, filtering for GPU/Subtimizer jobs."""
        try:
            # squeue format: %P (partition), %j (name)
            # Filter out interactive jobs (bash) or non-GPU jobs if possible
            cmd = f"squeue -u {self.user} -t RUNNING,PENDING -h --format='%P %j'"
            result = subprocess.run(
                cmd, 
                shell=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                text=True
            )
            
            lines = result.stdout.strip().split('\n')
            count = 0
            for line in lines:
                if not line.strip(): continue
                parts = line.split(maxsplit=1)
                if len(parts) < 2: continue
                
                partition, name = parts[0], parts[1]
                
                is_gpu = "gpu" in partition.lower()
                is_interactive = name == "bash"
                
                # Count if it's a GPU job AND not a generic bash session
                # Or if the name looks like a complex name (usually doesn't contain spaces)
                if is_gpu and not is_interactive:
                    count += 1
                elif "fold" in name.lower() or "design" in name.lower():
                     # Fallback if partition naming differs but job name is specific
                     count += 1
                     
            return count
        except Exception as e:
            print(f"Warning: Could not check job status ({e}). Assuming 0 jobs.")
            return 0
            
    def wait_for_slot(self, sleep_interval=15):
        """Blocks until there is a free slot (< max_jobs)."""
        while True:
            running = self.get_running_jobs()
            if running < self.max_jobs:
                return
            print(f"Queue full ({running}/{self.max_jobs}). Waiting {sleep_interval}s...")
            time.sleep(sleep_interval)
            
    def submit_job(self, script_path, job_name=None):
        """Submits a job via sbatch."""
        cmd = ["sbatch", script_path]
        if job_name:
            cmd.extend(["--job-name", job_name])
            
        print(f"Submitting: {script_path}")
        try:
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error submitting job: {e}")
            return False

def submit_slurm_job(job_name, script_content, work_dir):
    """
    Writes script_content to a file in work_dir and submits it via sbatch.
    Global utility to avoid instantiating JobManager for simple submissions.
    """
    script_path = os.path.join(work_dir, f"{job_name}.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
        
    print(f"Created submission script: {script_path}")
    
    cmd = ["sbatch", script_path]
    
    try:
        # Run sbatch from the work_dir
        subprocess.run(cmd, cwd=work_dir, check=True)
        print(f"Successfully submitted job: {job_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job {job_name}: {e}")
        return False
