import glob
import os
import shutil

def prepare_for_folding(file_path: str):
    """
    Prepares designed sequences for validation folding (af2-initial-guess)
    """
    print(f"Preparing sequences for folding from {file_path}")
    
    with open(file_path, 'r') as f:
        complexes = [line.strip() for line in f if line.strip()]

    for folder in complexes:
        print(f"Processing {folder}")
        
        # 1. Get Kinase Sequence from original fasta
        kinase_seq, header = _extract_kinase_seq(folder)
        if not kinase_seq:
            continue
            
        # 2. Setup paths
        af_complex = os.path.join(folder, "AFcomplex")
        clust_dir = os.path.join(af_complex, "mpnn_out_clust")
        fold_dir = os.path.join(af_complex, "mpnn_out_clust_fold")
        fold_seq_dir = os.path.join(fold_dir, "seqs")
        
        if os.path.exists(fold_dir):
            shutil.rmtree(fold_dir)
        os.makedirs(fold_seq_dir, exist_ok=True)
        
        # 3. Process clustered designs
        clustered_file = os.path.join(clust_dir, "all_design_clustered.fa")
        if not os.path.exists(clustered_file):
            print(f"  Missing clustered file: {clustered_file}")
            continue
            
        with open(clustered_file) as f:
            seq_count = 0
            for line in f:
                if ">" in line:
                    seq_count += 1
                    parts = line.strip().split(",")
                    des_id = f"des_{seq_count}"
                    
                    full_header = f">{seq_count}_{header}_{des_id}"
                    fname = full_header.replace(">", "")
                else:
                    peptide_seq = line.strip()
                    # Create directory for this specific design
                    des_seq_dir = os.path.join(fold_seq_dir, fname)
                    os.makedirs(des_seq_dir, exist_ok=True)
                    
                    # Write combined fasta
                    with open(os.path.join(des_seq_dir, f"{fname}.fasta"), "w") as out:
                        out.write(f"{full_header}\n")
                        out.write(f"{kinase_seq}:{peptide_seq}\n")
                        
        print(f"  Prepared {seq_count} sequences for folding.")

def _extract_kinase_seq(folder):
    """Finds the kinase sequence from the input fasta."""
    fastas = glob.glob(os.path.join(folder, "*.fasta"))
    if not fastas:
        print(f"  No fasta found in {folder}")
        return None, None
        
    with open(fastas[0]) as f:
        for line in f:
            if ">" in line:
                header = line.strip().replace(">", "")
            else:
                parts = line.split(":")
                return parts[0], header
    return None, None
