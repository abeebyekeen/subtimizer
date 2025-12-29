import os
import shutil

def setup_folders(file_path: str, setup_type: str):
    """
    Creates necessary folders based on the input file and setup type.
    
    Args:
        file_path: Path to the file containing complex/folder names.
        setup_type: Type of setup ('initial', 'mpnn', 'original').
    """
    
    # Read the list of folders/complexes
    try:
        with open(file_path, 'r') as f:
            folders = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    if setup_type == 'initial':
        _setup_initial(folders)
    elif setup_type == 'mpnn':
        _setup_mpnn(folders)
    elif setup_type == 'original':
        _setup_original(folders)

def _setup_initial(folders):
    count = 0
    for folder in folders:
        try:
            os.makedirs(folder, exist_ok=True)
            # Create subfolder required for AF
            os.makedirs(os.path.join(folder, 'AFcomplex'), exist_ok=True)
            count += 1
        except OSError as e:
            print(f"Error creating {folder}: {e}")
    print(f"Created {count} initial project folders.")

def _setup_mpnn(folders):
    count = 0 
    
    # Prepare JSON content
    import json
    config_data = {}
    
    for folder in folders:
        mpnn_path = os.path.join(folder, 'AFcomplex', 'mpnn_des')
        try:
            os.makedirs(mpnn_path, exist_ok=True)
            count += 1
            # Add default config for this complex
            config_data[folder] = {
                "chains_to_design": "B", 
                "fixed_positions": "4"
            }
        except OSError as e:
            print(f"Error creating MPNN folder for {folder}: {e}")
            
    # Write centralized config file if it doesn't exist
    config_file = "design_config.json"
    if not os.path.exists(config_file):
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=4)
        print(f"Created centralized configuration file: {config_file}")
    else:
        print(f"Note: {config_file} already exists. New complexes were NOT added to it.")
        
    print(f"Created {count} MPNN design folders.")

def _setup_original(folders):
    base_dir = "original_subs"
    os.makedirs(base_dir, exist_ok=True)
    
    count = 0
    for folder in folders:
        source_fastas = [f for f in os.listdir(folder) if f.endswith('.fasta')]
        
        target_dir = os.path.join(base_dir, folder)
        os.makedirs(target_dir, exist_ok=True)
        
        for fasta in source_fastas:
            shutil.copy(os.path.join(folder, fasta), os.path.join(target_dir, fasta))
            count += 1
            
    print(f"Setup {count} folders in {base_dir}.")
    
    print("Copying top ranked PDBs for original substrates...")
    pdb_count = 0
    import glob
    
    for folder in folders:
        # Source PDB: folder/AFcomplex/top5complex/*_relaxed_rank_001_*.pdb
        # Target: original_subs/folder/top5complex/
        
        # 1. Define Paths
        source_pdb_dir = os.path.join(folder, "AFcomplex", "top5complex")
        target_struct_dir = os.path.join(base_dir, folder, "top5complex")
        
        if not os.path.exists(source_pdb_dir):
            print(f"Warning: No structure found for {folder} (expected in {source_pdb_dir})")
            continue
            
        # 2. Find rank_001 pdb
        # Pattern: *_relaxed_rank_001_*.pdb
        pdbs = glob.glob(os.path.join(source_pdb_dir, "*_relaxed_rank_001_*.pdb"))
        if not pdbs:
            print(f"Warning: No rank_001 PDB found in {source_pdb_dir}")
            continue
            
        # 3. Create target dir
        os.makedirs(target_struct_dir, exist_ok=True)
        
        # 4. Copy
        for pdb in pdbs:
            try:
                shutil.copy(pdb, target_struct_dir)
                pdb_count += 1
            except Exception as e:
                print(f"Error copying PDB for {folder}: {e}")
                
    print(f"Copied {pdb_count} PDB structures for original validation to {base_dir}/<complex>/top5complex.")
