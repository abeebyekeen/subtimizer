
import os
import sys
import glob
import pandas as pd
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil



def submit_ipsae_job(complex_list_file, pae_cutoff, dist_cutoff, max_jobs, start=1, end=None):
    """Submits a SLURM job to run the ipSAE workflow."""
    from subtimizer.manager import submit_slurm_job
    from subtimizer.utils import get_template_content
    
    work_home = os.getcwd()
    abs_list_file = os.path.abspath(complex_list_file)
    
    template_content = get_template_content('ipsae_template.sh')

    # If end is None, pass -1 to the template to ensure a valid integer is passed to CLI
    end_val = end if end is not None else -1
        
    script_content = template_content.format(
        list_file=abs_list_file,
        pae_cutoff=pae_cutoff,
        dist_cutoff=dist_cutoff,
        max_jobs=max_jobs,
        start=start,
        end=end_val
    )
    
    job_name = f"ipsae_run_{start}_{end if end else 'end'}"
    submit_slurm_job(job_name, script_content, work_home)
    print(f"Submitted ipSAE job with {max_jobs} CPUs. Range: {start}-{end if end else 'end'}.")

def execute_ipsae_workflow(complex_list_file, pae_cutoff, dist_cutoff, max_jobs, start=1, end=None):
    """
    Main workflow for ipSAE evaluation (Internal Worker).
    Iterates through complexes, runs ipsae.py on PDBs/PAEs, and updates CSVs.
    """
    
    if end == -1:
        end = None
    
    ipsae_exe = shutil.which("ipsae.py")
    if not ipsae_exe:
        pass
        
    ipsae_exe = shutil.which("ipsae.py")
    if not ipsae_exe:
        common_paths = [
            os.path.join(os.getcwd(), "ipSAE", "ipsae.py"),
            os.path.join(os.environ.get("HOME", ""), "ipSAE", "ipsae.py")
        ]
        for p in common_paths:
            if os.path.exists(p):
                ipsae_exe = p
                break
    
    if not ipsae_exe:
        print("Error: 'ipsae.py' not found in PATH or standard locations.")
        sys.exit(1)

    work_home = os.getcwd()
    
    with open(complex_list_file, 'r') as f:
        all_complexes = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if end is None:
        complexes = all_complexes[start-1:]
    else:
        complexes = all_complexes[start-1:end]
        
    if not complexes:
        print(f"No complexes found in list range {start}-{end}.")
        return

    print(f"Processing {len(complexes)} complexes (Range: {start}-{end if end else 'end'})...")

    for complex_name in complexes:
        print(f"\nProcessing {complex_name}...")
        _process_complex(complex_name, work_home, ipsae_exe, pae_cutoff, dist_cutoff, max_jobs)

def _process_complex(complex_name, work_home, ipsae_exe, pae_cutoff, dist_cutoff, max_jobs):
    """Run ipSAE for a single complex."""
    csv_file = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
    
    if not os.path.exists(csv_file):
        print(f"  CSV not found for {complex_name}: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    tasks = []
    
    # Locate PDB and PAE files for each row
    for idx, row in df.iterrows():
        struct_id = str(row['id'])
        
        pae_file, pdb_file = _find_files(work_home, complex_name, struct_id, row)
        
        if pae_file and pdb_file:
            tasks.append((idx, pae_file, pdb_file))
        else:
            # Silence logging for skipped items
            pass

    if not tasks:
        print(f"  No valid PDB/PAE pairs found for {complex_name}.")
        return

    print(f"  Running ipSAE on {len(tasks)} structures using {max_jobs} workers...")
    
    results = {}
    with ProcessPoolExecutor(max_workers=int(max_jobs)) as executor:
        futures = {executor.submit(_run_single_ipsae, t, ipsae_exe, pae_cutoff, dist_cutoff): t[0] for t in tasks}
        
        for future in as_completed(futures):
            idx, status, output = future.result()
            if status == "SUCCESS":
                results[idx] = output # output is pdb_path (used to find txt output)
            else:
                print(f"    Row {idx} Failed: {output}")

    # Parse outputs and update CSV
    _update_csv_with_results(df, results, pae_cutoff, dist_cutoff, csv_file)

def _find_files(work_home, complex_name, struct_id, row):
    """Locate PDB and PAE files."""
    
    base_dirs = [
        # Standard Designed
        os.path.join(work_home, complex_name, "AFcomplex", "mpnn_out_clust_fold", "seqs", struct_id),
        # Flat original
        os.path.join(work_home, complex_name, "top5complex"), 
    ]
    
    pae_file = None
    pdb_file = None
    
    # Try Standard
    search_dir = base_dirs[0]
    if os.path.exists(search_dir):
        # find JSON
        jsons = glob.glob(os.path.join(search_dir, f"{struct_id}*_rank_001_*.json"))
        pdbs = glob.glob(os.path.join(search_dir, f"{struct_id}*rank_001_*.pdb"))
        if jsons and pdbs:
            return jsons[0], pdbs[0]

    # Try Flat Original (top5complex)
    flat_dir = base_dirs[1]
    if os.path.exists(flat_dir):
        pass

    # Fallback: Validation Round folders (Legacy special case)
    # Check AFcomplex/round_*/
    for r in range(1, 6):
        r_dir = os.path.join(work_home, complex_name, "AFcomplex", f"round_{r}")
        if os.path.exists(r_dir):
            # Look for file matching ID
            fold_val = str(row.get('fold', struct_id)) 
            
            jsons = glob.glob(os.path.join(r_dir, f"*{fold_val}*.json"))
            pdbs = glob.glob(os.path.join(r_dir, f"*{fold_val}*.pdb"))
            
            if jsons and pdbs:
                 # Prioritize relaxed
                 relaxed = [p for p in pdbs if "relaxed" in p and "unrelaxed" not in p]
                 final_pdb = relaxed[0] if relaxed else pdbs[0]
                 return jsons[0], final_pdb

    return None, None

def _run_single_ipsae(task, ipsae_exe, pae_cutoff, dist_cutoff):
    idx, pae, pdb = task
    
    # Command: python ipsae.py <pae> <pdb> <pae_cut> <dist_cut>
    # Note: ipsae.py generates output in same dir as pdb
    cmd = ["python3" if ipsae_exe.endswith(".py") else ipsae_exe, pae, pdb, str(pae_cutoff), str(dist_cutoff)]
    if ipsae_exe.endswith(".py"):
        cmd = ["python3", ipsae_exe, pae, pdb, str(pae_cutoff), str(dist_cutoff)]
    else:
        cmd = [ipsae_exe, pae, pdb, str(pae_cutoff), str(dist_cutoff)]

    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            return idx, "SUCCESS", pdb
        else:
            return idx, "ERROR", res.stderr
    except Exception as e:
        return idx, "EXCEPTION", str(e)

def _update_csv_with_results(df, results, pae_cutoff, dist_cutoff, csv_file):
    
    new_cols = ['ipSAE', 'ipSAE_d0chn', 'ipSAE_d0dom', 'ipTM_af', 'ipTM_d0chn', 'ipSAE_min']
    for c in new_cols:
        if c not in df.columns: df[c] = None # init
        
    p_str = str(int(pae_cutoff)).zfill(2)
    d_str = str(int(dist_cutoff)).zfill(2)
    
    update_count = 0
    for idx, pdb_path in results.items():
        base = os.path.splitext(pdb_path)[0]
        txt_path = f"{base}_{p_str}_{d_str}.txt"
        
        if os.path.exists(txt_path):
            parsed = _parse_ipsae_txt(txt_path)
            if parsed:
                for k, v in parsed.items():
                    if k in new_cols:
                        df.at[idx, k] = v
                update_count += 1
                
    out_file = csv_file.replace(".csv", f"_with_ipSAEmin_{pae_cutoff}_{dist_cutoff}.csv")
    df.to_csv(out_file, index=False)
    print(f"  Updated CSV saved to: {out_file}")
    print(f"  Updated {update_count} entries.")

def _parse_ipsae_txt(txt_path):
    data = {}
    ipsae_vals = []
    try:
        with open(txt_path, 'r') as f:
            for line in f:
                parts = line.split()
                if not parts: continue
                
                # Check metrics
                if "max" in parts:
                    try:
                        i = parts.index("max")
                        # parts: ... max ipSAE ipSAE_d0chn ...
                        data['ipSAE'] = float(parts[i+1])
                        data['ipSAE_d0chn'] = float(parts[i+2])
                        data['ipSAE_d0dom'] = float(parts[i+3])
                        data['ipTM_af'] = float(parts[i+4])
                        data['ipTM_d0chn'] = float(parts[i+5])
                        ipsae_vals.append(data['ipSAE'])
                    except: pass
                elif "asym" in parts:
                    try:
                        i = parts.index("asym")
                        ipsae_vals.append(float(parts[i+1]))
                    except: pass
        
        if ipsae_vals:
            data['ipSAE_min'] = min(ipsae_vals)
            
        return data if 'ipSAE' in data else None
    except:
        return None
