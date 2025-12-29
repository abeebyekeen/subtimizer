import os
import glob
from Bio import PDB
from subtimizer.manager import JobManager

import pkgutil

def run_pdb_fix(file_path: str, max_jobs: int = 4, start: int = 1, end: int = None):
    """
    PDB fixing for AF2 initial guess
    """
    print(f"Running PDB fixing for complexes in {file_path}")
    
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
            init_guess_in = os.path.join(work_dir, "af2_init_guess_in")
            structs_dir = os.path.join(work_dir, "structs")
        elif os.path.exists(flat_top5_dir):
            mode = "flat"
            work_dir = complex_name
            init_guess_in = os.path.join(work_dir, "af2_init_guess_in")
            structs_dir = flat_top5_dir
        else:
            print(f"Skipping {complex_name}: No valid structure directory found (checked Standard and Flat paths).")
            continue

        # Clean and recreate input dir
        if os.path.exists(init_guess_in):
            for f in glob.glob(os.path.join(init_guess_in, "*")):
                try:
                    os.remove(f)
                except:
                    pass
        else:
            os.makedirs(init_guess_in, exist_ok=True)
            
        # Search for PDBs
        # Mode determines primary source, but we use the set path
        pdbs = glob.glob(os.path.join(structs_dir, "*_rank_001_*.pdb"))
        if not pdbs:
             pdbs = glob.glob(os.path.join(structs_dir, "*.pdb"))
             
        if not pdbs:
            print(f"No PDBs found in {structs_dir} ({mode} mode).")
            continue
            
        # Copy PDBs
        import shutil
        for pdb in pdbs:
            shutil.copy(pdb, init_guess_in)

        # Create worker script
        script_path = os.path.join(work_dir, "run_pdb_fix.sh")
        _write_fix_script(script_path, complex_name, os.path.abspath(init_guess_in))
        
        wd = os.getcwd()
        os.chdir(work_dir)
        manager.submit_job("run_pdb_fix.sh", job_name=f"fix_{complex_name}")
        os.chdir(wd)

def _write_fix_script(path, complex_name, target_dir):
    """
    Writes a python script wrapper using the external template.
    """
    try:
        template_bytes = pkgutil.get_data('subtimizer.templates', 'fix_pdb_template.sh')
        if template_bytes is None:
             # Fallback if template missing (though it shouldn't be)
             raise FileNotFoundError("fix_pdb_template.sh not found")
        template_content = template_bytes.decode('utf-8')
    except Exception as e:
        print(f"Error loading template: {e}")
        return

    content = template_content.format(
        complex_name=complex_name,
        target_dir=target_dir
    )
    
    with open(path, 'w') as f:
        f.write(content)

def fix_pdbs_in_dir(directory: str):
    """
    actual logic to fix PDBs using BioPython
    """
    print(f"Processing PDBs in {directory}")
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()
    
    for pdb_file in glob.glob(os.path.join(directory, "*.pdb")):
        try:
            structure = parser.get_structure('struct', pdb_file)
            
            # Logic: Ensure Chain IDs are standard (A, B) and Residues are renumbered?
            # performing a standard renumbering:
            # Chain 1 -> A, Chain 2 -> B
            # Residues 1..N
            
            new_chain_ids = ['A', 'B', 'C', 'D', 'E']
            # Logic: Ensure Binder (shorter chain) is Chain A, Kinase (longer) is Chain B.
            # and no overlapping residues.
            
            chains = list(structure.get_chains())
            
            # Sort chains by length (number of residues)
            # Shortest first -> Chain A (Peptide), Longest -> Chain B (Kinase)
            chains.sort(key=lambda c: len(list(c)), reverse=False)
            
            # Helper to get length safely (len() works on Chain)
            def get_len(c): return len(list(c))

            kinase_len = 0
            if len(chains) > 1:
                # Assuming second chain is Kinase (longer)
                kinase_len = get_len(chains[1])
                        
            start_nums = []
            
            peptide_start = 801
            if kinase_len >= 800:
                peptide_start = kinase_len + 1
                
            start_nums.append(peptide_start) # For Chain A (Peptide)
            start_nums.append(1)             # For Chain B (Kinase)
            
            current_max = max(peptide_start + get_len(chains[0]), 1 + kinase_len)

            for k in range(2, len(chains)):
                start_nums.append(current_max + 1)
                current_max += get_len(chains[k])
                
            # Determine final mappings (Old Chain Object -> New ID, New Start Res)
            chain_updates = []
            
            peptide_start = 801
            if kinase_len >= 800:
                peptide_start = kinase_len + 1
            
            import random
            
            # Assign Temp IDs
            for i, chain in enumerate(chains):
                chain.id = f"T{i}" 
            
            # Second pass: Assign Final IDs and Renumber
            current_max = max(peptide_start + get_len(chains[0]), 1 + kinase_len)
            
            for i, chain in enumerate(chains):
                # Target ID
                new_id = new_chain_ids[i] if i < len(new_chain_ids) else f"X{i}"
                
                # Target Start
                if i == 0: start_res = peptide_start
                elif i == 1: start_res = 1
                else: 
                     start_res = current_max + 1
                     current_max += get_len(chain)
                
                chain.id = new_id
                
                for j, residue in enumerate(chain):
                     residue.id = (' ', start_res + j, ' ')
                     
            io.set_structure(structure)
            io.save(pdb_file)
            # print(f"Fixed {os.path.basename(pdb_file)}")
            
        except Exception as e:
            print(f"Error fixing {pdb_file}: {e}")
