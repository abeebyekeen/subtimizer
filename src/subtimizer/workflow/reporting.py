import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import colorcet as cc
from scipy import stats

def run_reporting(file_path: str, start: int = 1, end: int = None):
    """
    data merging and advanced plotting
    """
    print(f"Generating Reports for complexes in {file_path}")
    
    with open(file_path, 'r') as f:
        all_complexes = [line.strip() for line in f if line.strip()]

    # Filter by start/end
    if end is None:
        complexes = all_complexes[start-1:]
    else:
        complexes = all_complexes[start-1:end]
        
    if not complexes:
        print(f"No complexes found in range {start}-{end if end else 'end'}.")
        return

    print(f"Processing {len(complexes)} complexes (Range: {start}-{end if end else 'end'})...")

    _merge_metrics(complexes)
    
    _merge_with_original(complexes)

    _add_sequences(complexes)
    
    _plot_scatter(complexes, file_path, start, end)

    _plot_swarm(complexes, file_path, start, end)

    _plot_ipsae(complexes, file_path, start, end)

    _plot_scatter_ipsae_colored(complexes, file_path, start, end)

def _merge_metrics(complexes):
    """
    Merges AF2 scores (ipTM, pTM) with folding logs.
    """
    print("MetaData Merging...")
    work_home = os.getcwd()
    
    for complex_name in complexes:
        # Paths
        af_complex_dir = os.path.join(complex_name, "AFcomplex")
        fold_dir = os.path.join(af_complex_dir, "mpnn_out_clust_fold")
        # Try Standard Path
        init_guess_dir = os.path.join(fold_dir, "af2_init_guess.rec8")
        is_flat = False
        
        if not os.path.exists(init_guess_dir):
            flat_guess = os.path.join(complex_name, "af2_init_guess.rec8")
            if os.path.exists(flat_guess):
                init_guess_dir = flat_guess
                is_flat = True
        
        score_file = os.path.join(init_guess_dir, "af2score.dat")
        if not os.path.exists(score_file):
            print(f"Skipping {complex_name}: af2score.dat not found.")
            continue
            
        merged_data = []
        try:
            with open(score_file) as sf:
                lines = sf.readlines()
                
            for i, line in enumerate(lines):
                parts = line.strip().split()
                if i == 0:
                   # Header
                   headers = ["id", "fold", "ipTM", "pTM", "pLDDT"] + parts[1:-1] # Reordered to match requirement
                   merged_data.append(headers)
                else:
                    # Data Line: SCORE: <metrics> description
                    description = parts[-1] 
                    # Desc format example: "design_001_unrelaxed_..."
                    design_id = description.split("_unrelaxed_")[0]
                    
                    # fetch folding log for this design
                    log_file = os.path.join(fold_dir, "seqs", design_id, "log.txt")
                    
                    found_log = False
                    if os.path.exists(log_file):
                         found_log = True
                    else:
                        parent_project_dir = os.path.abspath(os.path.join(work_home, "..", complex_name))
                        
                        search_roots = []
                        if "original_subs" in work_home:
                             search_roots.append(os.path.join(work_home, "..", complex_name))
                             
                        for root in search_roots:
                            # Search rounds 1-5
                            for r in range(1, 6):
                                possible_log = os.path.join(root, "AFcomplex", f"round_{r}", "log.txt")
                                if os.path.exists(possible_log):
                                    log_file = possible_log
                                    found_log = True
                                    break
                            if found_log: break
                    
                    # Default values if log missing
                    iptm, ptm, plddt, fld = "N/A", "N/A", "N/A", "N/A"
                    
                    if found_log and os.path.exists(log_file):
                        with open(log_file) as lf:
                            for logline in lf:
                                # For original substrates, strictly look for rank_001
                                if "rank_001" in logline:
                                    lp = logline.strip().split()
                                    try:
                                        iptm = lp[-1].replace("ipTM=", "")
                                        ptm = lp[-2].replace("pTM=", "")
                                        plddt = lp[-3].replace("pLDDT=", "")
                                        fld = lp[-4].split("_v3_")[0] # model_1_multimer
                                    except:
                                        pass # Keep N/A
                                    break
                    
                    # Combine
                    # New Order: id, fold, ipTM, pTM, pLDDT, ...
                    row = [design_id, fld, iptm, ptm, plddt] + parts[1:-1] 
                    merged_data.append(row)
            
            import pandas as pd
            df = pd.DataFrame(merged_data[1:], columns=merged_data[0])
            
            if "pTM" in df.columns and "ipTM" in df.columns:
                try:
                    df["pTM"] = pd.to_numeric(df["pTM"], errors='coerce')
                    df["ipTM"] = pd.to_numeric(df["ipTM"], errors='coerce')
                    
                    # 0.2 * pTM + 0.8 * ipTM
                    df["pTM_ipTM"] = (0.2 * df["pTM"] + 0.8 * df["ipTM"]).round(3)
                    
                    # Move to position 5 if possible 
                    cols = list(df.columns)
                    cols.insert(5, cols.pop(cols.index("pTM_ipTM")))
                    df = df[cols]
                except Exception as e:
                    print(f"Warning: Could not calculate pTM_ipTM for {complex_name}: {e}")

            out_csv = os.path.join(init_guess_dir, f"{complex_name}_merged_scores_pTM-ipTM.csv")
            df.to_csv(out_csv, index=False)
            
            central_dest = os.path.join(work_home, "af2_init_guess", "data", complex_name)
            os.makedirs(central_dest, exist_ok=True)
            try:
                import shutil
                shutil.copy2(out_csv, central_dest)
            except Exception as e:
                 print(f"Warning: Could not copy report to {central_dest}: {e}")
            
        except Exception as e:
            print(f"Error merging {complex_name}: {e}")

def _merge_with_original(complexes):
    """
    Merges designed substrate data with original substrate data (if available).
    """
    if "original_subs" in os.getcwd():
        return

    orig_subs_path = "original_subs"
    if not os.path.exists(orig_subs_path):
        return

    print("Merging with Original Substrates data...")
    work_home = os.getcwd()

    for complex_name in complexes:
        designed_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM.csv")
        
        orig_csv_source = os.path.join(orig_subs_path, complex_name, "AFcomplex", "mpnn_out_clust_fold", 
                                       "af2_init_guess.rec8", f"{complex_name}_merged_scores_pTM-ipTM.csv")
                                       
        if not os.path.exists(orig_csv_source):
            orig_csv_source = os.path.join(orig_subs_path, complex_name, "af2_init_guess.rec8", 
                                           f"{complex_name}_merged_scores_pTM-ipTM.csv")
                                       
        if not os.path.exists(designed_csv):
            # No designed data? Skip
            continue
            
        if not os.path.exists(orig_csv_source):
            print(f"Info: Original data not found for {complex_name} (checked {orig_csv_source}). Skipping merge.")
            continue
            
        try:
            df_des = pd.read_csv(designed_csv)
            df_orig = pd.read_csv(orig_csv_source)
            
            if not df_orig.empty:
                # Take top row only
                df_orig_top = df_orig.iloc[[0]]
                df_final = pd.concat([df_des, df_orig_top], ignore_index=True)
            else:
                df_final = df_des
            
            out_file = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
            df_final.to_csv(out_file, index=False)
            print(f"Created combined report: {out_file}")
            
        except Exception as e:
            print(f"Error merging original data for {complex_name}: {e}")

def _plot_swarm(complexes, source_file, start=1, end=None):
    """
    Generates Swarm Plots
    """
    print("Generating Swarm Plots...")

    plt.rcParams['axes.linewidth'] = 1.7
    plt.rcParams['xtick.major.width'] = 1.7
    plt.rcParams['ytick.major.width'] = 1.7
    plt.rcParams['xtick.major.size'] = 5.2
    plt.rcParams['ytick.major.size'] = 4.0

    palette = sns.color_palette(cc.glasbey, n_colors=len(complexes)+10)

    plt.figure(figsize=(6, 10))
    
    plotted_count = 0
    
    for idx, complex_name in enumerate(complexes):
        work_home = os.getcwd()
        combined_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
        
        fallback_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM.csv")
        
        csv_file = combined_csv if os.path.exists(combined_csv) else fallback_csv
        
        if not os.path.exists(csv_file):
            csv_file = os.path.join(work_home, complex_name, "AFcomplex", "mpnn_out_clust_fold", 
                                  "af2_init_guess.rec8", f"{complex_name}_merged_scores_pTM-ipTM.csv")
            if not os.path.exists(csv_file):
                 csv_file = os.path.join(work_home, complex_name, "af2_init_guess.rec8", 
                                        f"{complex_name}_merged_scores_pTM-ipTM.csv")
                 if not os.path.exists(csv_file):
                    continue
        
        # Check if this is the combined file (implies parental data exists)
        is_combined = "with_oriSubs.csv" in csv_file

        try:
            df = pd.read_csv(csv_file)
            
            metric_col = None
            if 'pae_interaction' in df.columns:
                metric_col = 'pae_interaction'
            elif 'ipAE' in df.columns:
                metric_col = 'ipAE'
            else:
                # Fallback to column 7 if it exists
                if len(df.columns) > 7:
                    metric_col = df.columns[7]
            
            if not metric_col:
                print(f"Warning: Could not identify plotting column for {complex_name}")
                continue
                
            data_values = df[metric_col]
            
            # Ensure numeric
            df[metric_col] = pd.to_numeric(df[metric_col], errors='coerce')
            df = df.dropna(subset=[metric_col])
            
            data_values = df[metric_col]
            
            if len(data_values) == 0:
                continue

            plotted_count += 1
            
            p_name = complex_name
            
            # Choose Swarm Color (Avoid Red-ish)
            # if RGB has high Red component, pick another.
            swarm_color = palette[idx % len(palette)]
            if len(swarm_color) == 3 and swarm_color[0] > 0.7 and swarm_color[1] < 0.5:
                 # Shift by 10 or something arbitrary to get away from red
                 swarm_color = palette[(idx + 10) % len(palette)]
                            
            if is_combined:
                # Plot Design (Swarm)
                # Assumes Last Row is Parental
                if len(data_values) > 1:
                    sns.swarmplot(y=[p_name] * (len(data_values) - 1), x=data_values[:-1], 
                                 alpha=0.9, color=swarm_color, size=6)
                
                # Plot Parental (Last Point) - Red Dot
                plt.scatter(y=p_name, x=data_values.iloc[-1], s=60, color='red', edgecolor='black', zorder=28)
            else:
                # Plot All as Swarm (No Red Dot)
                sns.swarmplot(y=[p_name] * len(data_values), x=data_values, 
                             alpha=0.9, color=swarm_color, size=6)

        except Exception as e:
            print(f"Error plotting {complex_name}: {e}")

    if plotted_count == 0:
        print("No data found for plotting.")
        return

    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.ylabel('Kinase-peptide complex', fontsize=20, fontweight='bold')
    plt.xlabel('ipAE' + r' ($\AA$)', fontsize=20, fontweight='bold')
    
    ax = plt.gca()
    ax.tick_params(axis='x', which='major', length=5.2, width=1.7)
    ax.tick_params(axis='x', which='minor', length=4, width=1.4)
    
    base_name = os.path.basename(source_file)
    file_stem = os.path.splitext(base_name)[0]
    end_val = end if end else "end"
    out_file = f"validation_swarmplot_{file_stem}_{start}_{end_val}.png"
    
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {out_file}")
    plt.clf()

def _add_sequences(complexes):
    """
    Adds peptide sequence column to the final CSVs.
    """
    print("Adding sequences to reports...")
    from Bio import SeqIO
    # Helper to get seq from PDB
    def get_seq_from_pdb(pdb_path):
        try:
            from Bio.PDB import PDBParser
            from Bio.SeqUtils import seq1
            
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('struct', pdb_path)
            model = structure[0]
            
            shortest_seq = None
            min_len = float('inf')
            
            for chain in model:
                seq = ""
                for res in chain:
                    if res.id[0] == " ": 
                        try:
                            rname = res.get_resname()
                            s = seq1(rname)
                            if s: seq += s
                        except:
                            pass
                
                if len(seq) > 0:
                    if len(seq) < min_len:
                        min_len = len(seq)
                        shortest_seq = seq
            
            if not shortest_seq:
                print(f"Debug: No protein chain found in {pdb_path}")
                return "N/A"
            return shortest_seq
            
        except Exception as e:
             print(f"Debug: Error parsing {pdb_path}: {e}")
             return "N/A"

    work_home = os.getcwd()

    for complex_name in complexes:
        # Paths setup
        combined_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
        fallback_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM.csv")
        csv_file = combined_csv if os.path.exists(combined_csv) else fallback_csv
        
        if not os.path.exists(csv_file):
             continue

        try:
             df = pd.read_csv(csv_file)
             
             seqs = []
             for i, row in df.iterrows():
                 d_id = str(row.get('id', ''))
                 
                 # 1. Try finding Folded PDB (Standard)
                 # Try relaxed first, then unrelaxed
                 patterns = [
                     f"{d_id}_relaxed_rank_001_*.pdb",
                     f"{d_id}_unrelaxed_rank_001_*.pdb",
                     f"{d_id}_rank_001_*.pdb"
                 ]
                 found_seq = False
                 for pat in patterns:
                     pdb_path = os.path.join(work_home, complex_name, "AFcomplex", "mpnn_out_clust_fold", "structs", pat)
                     matches = glob.glob(pdb_path)
                     if matches:
                         seqs.append(get_seq_from_pdb(matches[0]))
                         found_seq = True
                         break
                 if found_seq: continue
                 
                 # 2. Flat Mode
                 matches_flat = glob.glob(os.path.join(work_home, complex_name, "top5complex", "*.pdb"))
                 if matches_flat and (d_id in complex_name or complex_name in d_id):
                      seqs.append(get_seq_from_pdb(matches_flat[0]))
                      continue
                 
                 # 3. Original Substrate Row (Last Row Check)
                 if i == len(df) - 1:
                     # finding any PDB in standard Original path
                     std_orig = os.path.join(work_home, "original_subs", complex_name, "AFcomplex", "top5complex", "*.pdb")
                     matches_std = glob.glob(std_orig)
                     if matches_std:
                         seqs.append(get_seq_from_pdb(matches_std[0]))
                         continue
                     
                     # Try finding any PDB in Flat Mode (if running from main dir looking into original_subs)
                     flat_orig = os.path.join(work_home, "original_subs", complex_name, "top5complex", "*.pdb")
                     matches_flat_orig = glob.glob(flat_orig)
                     if matches_flat_orig:
                         seqs.append(get_seq_from_pdb(matches_flat_orig[0]))
                         continue

                 # Debug failure
                 print(f"Debug: Could not find PDB for ID {d_id} in {complex_name}")
                 seqs.append("N/A")
             
             df['pep_sequence'] = seqs
             df.to_csv(csv_file, index=False)
             print(f"Updated {csv_file} with sequences.")

        except Exception as e:
             print(f"Error adding sequences for {complex_name}: {e}")

def _plot_scatter(complexes, source_file, start=1, end=None):
    """
    Generates Scatter Plots matching legacy Script 33.
    """
    print("Generating Scatter Plots...")
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.linewidth'] = 1.5
    
    LABEL_SIZE = 21
    TICK_SIZE = 20
    LEGEND_SIZE = 14
    
    for complex_name in complexes:
         work_home = os.getcwd()
         combined_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
         fallback_csv = os.path.join(work_home, "af2_init_guess", "data", complex_name, f"{complex_name}_merged_scores_pTM-ipTM.csv")
         csv_file = combined_csv if os.path.exists(combined_csv) else fallback_csv
         
         if not os.path.exists(csv_file):
             continue
             
         output_dir = os.path.dirname(csv_file)
             
         try:
             df = pd.read_csv(csv_file)
             
             # Metric columns
             metrics = {}
             # Find ipAE
             if 'ipAE' in df.columns: metrics['y'] = 'ipAE'
             elif 'pae_interaction' in df.columns: metrics['y'] = 'pae_interaction'
             else: metrics['y'] = df.columns[7] if len(df.columns) > 7 else None
             
             # Find pTM_ipTM
             if 'pTM_ipTM' in df.columns: metrics['x'] = 'pTM_ipTM'
             else: metrics['x'] = None 
             
             if not metrics['x'] or not metrics['y']:
                 print(f"Skipping scatter for {complex_name}: missing columns.")
                 continue
                 
             # Verify dtypes
             df[metrics['x']] = pd.to_numeric(df[metrics['x']], errors='coerce')
             df[metrics['y']] = pd.to_numeric(df[metrics['y']], errors='coerce')
             df = df.dropna(subset=[metrics['x'], metrics['y']])
             
             if df.empty: continue
             
             plddt_col = None
             for c in ['plddt_binder', 'pLDDT_binder', 'pLDDT']:
                 if c in df.columns:
                     plddt_col = c
                     break
             if plddt_col:
                 df[plddt_col] = pd.to_numeric(df[plddt_col], errors='coerce')

             def setup_axis(ax):
                 ax.set_ylabel(r'$\mathbf{ipAE}$ ($\AA$)', fontsize=LABEL_SIZE)
                 ax.set_xlabel(r'$\mathbf{0.2pTM+0.8ipTM}$', fontsize=LABEL_SIZE)
                 ax.tick_params(axis='both', which='major', labelsize=TICK_SIZE, length=6, width=1.7)
                 ax.tick_params(axis='x', which='minor', length=4, width=1.5)
                 
                 for label in ax.get_xticklabels() + ax.get_yticklabels():
                     label.set_fontweight('bold')
                     
                 ax.yaxis.set_major_locator(ticker.MultipleLocator(6))
                 ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
                 ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
                 
                 ax.set_xlim(0.35, 0.96)
                 ax.set_ylim(5.7, 30)
                 
                 for spine in ax.spines.values():
                    spine.set_linewidth(1.5)
                 
                 ax.set_title(f'{complex_name}', fontsize=10)

             is_combined = "with_oriSubs" in csv_file
             
             size = 100
             size_star = 600
             star_edge = 4
             
             # --- Plot 1: Simple Scatter ---
             fig, ax = plt.subplots(figsize=(6, 6))
             
             if is_combined:
                 ax.scatter(x=df[metrics['x']][:-1], y=df[metrics['y']][:-1], s=size)
                 ax.scatter(x=df[metrics['x']].iloc[-1], y=df[metrics['y']].iloc[-1], 
                            marker='*', s=size_star, c='red', edgecolor='black', linewidths=star_edge, zorder=10, label='Parental')
                 ax.legend(loc='lower left', fontsize=LEGEND_SIZE)
             else:
                 ax.scatter(x=df[metrics['x']], y=df[metrics['y']], s=size)
            
             setup_axis(ax)
             figname = os.path.join(output_dir, f"{complex_name}_scatter_pTM-ipTM_vs_ipAE.png")
             plt.savefig(figname, dpi=300, bbox_inches='tight')
             plt.close()
             
             # --- Plot 2: Colored by pLDDT ---
             if plddt_col:
                 fig, ax = plt.subplots(figsize=(7, 6))
                 cmap = plt.get_cmap('turbo_r')
                 norm = plt.Normalize(vmin=28, vmax=80) 
                 
                 if is_combined:
                     # Designs
                     ax.scatter(x=df[metrics['x']][:-1], y=df[metrics['y']][:-1], 
                                     c=df[plddt_col][:-1], cmap=cmap, norm=norm, s=size)
                     # Parental
                     ax.scatter(x=df[metrics['x']].iloc[-1], y=df[metrics['y']].iloc[-1], 
                                marker='*', s=size_star, c=[df[plddt_col].iloc[-1]], cmap=cmap, norm=norm,
                                edgecolor='black', linewidths=star_edge, zorder=10, label='Parental')
                     ax.legend(loc='lower left', fontsize=LEGEND_SIZE)
                 else:
                     ax.scatter(x=df[metrics['x']], y=df[metrics['y']], 
                                     c=df[plddt_col], cmap=cmap, norm=norm, s=size)
                     
                 # Colorbar
                 sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                 cbar = plt.colorbar(sm, ax=ax)
                 cbar.set_label(r'$\mathbf{pLDDT}$' + r'_$\mathbf{peptide}$', fontsize=LABEL_SIZE)
                 cbar.ax.tick_params(labelsize=19)
                 for l in cbar.ax.get_yticklabels(): l.set_fontweight('bold')
                 
                 setup_axis(ax)
                 figname2 = os.path.join(output_dir, f"{complex_name}_scatter_pTM-ipTM_vs_ipAE_vs_pLDDT.png")
                 plt.savefig(figname2, dpi=300, bbox_inches='tight')
                 plt.close()
                 
                 # --- Plot 3: Labelled (Top 5 lowest ipAE) ---
                 fig, ax = plt.subplots(figsize=(7, 6))
                 
                 # Same scatter
                 if is_combined:
                     ax.scatter(x=df[metrics['x']][:-1], y=df[metrics['y']][:-1], 
                                c=df[plddt_col][:-1], cmap=cmap, norm=norm, s=size)
                     ax.scatter(x=df[metrics['x']].iloc[-1], y=df[metrics['y']].iloc[-1], 
                                marker='*', s=size_star, c=[df[plddt_col].iloc[-1]], cmap=cmap, norm=norm,
                                edgecolor='black', linewidths=star_edge, zorder=10, label='Parental')
                     ax.legend(loc='lower left', fontsize=LEGEND_SIZE)
                 else:
                     ax.scatter(x=df[metrics['x']], y=df[metrics['y']], 
                                c=df[plddt_col], cmap=cmap, norm=norm, s=size)
                                
                 cbar = plt.colorbar(sm, ax=ax)
                 cbar.set_label(r'$\mathbf{pLDDT}$' + r'_$\mathbf{peptide}$', fontsize=LABEL_SIZE)
                 cbar.ax.tick_params(labelsize=19)
                 for l in cbar.ax.get_yticklabels(): l.set_fontweight('bold')
                 
                 setup_axis(ax)
                 
                 # Labelling Top 5 Lowest ipAE
                 lowest = df.nsmallest(5, metrics['y'])
                 for _, row in lowest.iterrows():
                    rid = str(row['id'])
                    label = rid
                    if "design" in rid:
                        idx = rid.find("design")
                        parts = rid[idx:].split("_")
                        if len(parts) > 1: label = f"d{parts[1]}" # d001
                     
                    plt.text(row[metrics['x']], row[metrics['y']], label, fontsize=10, fontweight='bold')
                 
                 figname3 = os.path.join(output_dir, f"{complex_name}_scatter_pTM-ipTM_vs_ipAE_vs_pLDDT_labelled.png")
                 plt.savefig(figname3, dpi=300, bbox_inches='tight')
                 plt.close()
                 
                 print(f"Generated scatter plots for {complex_name} in {output_dir}")
                 
         except Exception as e:
             print(f"Error scatter plotting {complex_name}: {e}")

def _plot_ipsae(complexes, source_file, start=1, end=None):
    """
    Generates ipSAE correlation plots
    """
    print("Checking for ipSAE data to plot...")
    
    dot_size = 105
    fig_size = (16, 18)
    TITLE_SIZE = 22
    LABEL_SIZE = 23
    TEXT_SIZE = 24
    FONT_SIZE = 22
    
    def get_stats(x, y):
        if len(x) < 2: return np.nan, np.nan, np.nan, np.nan
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        return slope, intercept, r_value**2, r_value

    def apply_custom_style(ax, x_ticks_cfg=None, y_ticks_cfg=None):
        if x_ticks_cfg:
            ax.set_xticks(x_ticks_cfg['major'])
            ax.set_xticklabels([f"{t:.1f}" if isinstance(t, float) else str(int(t)) for t in x_ticks_cfg['major']], 
                               fontsize=FONT_SIZE, fontweight='bold')
            if x_ticks_cfg.get('minor') is not None: ax.set_xticks(x_ticks_cfg['minor'], minor=True)
                
        if y_ticks_cfg:
            ax.set_yticks(y_ticks_cfg['major'])
            ax.set_yticklabels([f"{t:.1f}" if isinstance(t, float) else str(int(t)) for t in y_ticks_cfg['major']], 
                               fontsize=FONT_SIZE, fontweight='bold')
            if y_ticks_cfg.get('minor') is not None: ax.set_yticks(y_ticks_cfg['minor'], minor=True)

        ax.tick_params(axis='both', which='major', length=7, width=1.8, labelsize=FONT_SIZE)
        ax.tick_params(axis='x', which='minor', length=5.5, width=1.6)
        ax.tick_params(axis='y', which='minor', length=5.5, width=1.6)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(FONT_SIZE)
            label.set_fontweight('bold')
        for spine in ax.spines.values(): spine.set_linewidth(1.7)

    for complex_name in complexes:
        work_home = os.getcwd()
        data_dir = os.path.join(work_home, "af2_init_guess", "data", complex_name)
        patterns = [
            os.path.join(data_dir, f"*_with_ipSAEmin_*.csv"),
            os.path.join(data_dir, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
        ]
         
        target_csv = None
        for p in patterns:
            matches = glob.glob(p)
            if matches:
                # Check columns
                for m in matches:
                    try:
                        d = pd.read_csv(m)
                        if 'ipSAE' in d.columns or 'ipSAE_min' in d.columns:
                            target_csv = m
                            break
                    except: pass
            if target_csv: break
             
        if not target_csv:
            continue
             
        print(f"Generating ipSAE plots for {complex_name}...")
        df = pd.read_csv(target_csv)
        output_dir = os.path.dirname(target_csv)
        base_name = os.path.splitext(os.path.basename(target_csv))[0]
         
        mapper = {'pae_interaction': 'ipAE', 'plddt_binder': 'pLDDT_peptide'}
        df = df.rename(columns={k:v for k,v in mapper.items() if k in df.columns})
         
        # Define plotting tasks: (y_metric, filename_suffix)
        tasks = []
        if 'ipSAE' in df.columns: tasks.append(('ipSAE', 'ipSAE'))
        if 'ipSAE_min' in df.columns: tasks.append(('ipSAE_min', 'ipSAEmin'))
         
        for y_metric, suffix in tasks:
            # Configuration
            pairs_config = [
                ('ipTM', y_metric),           
                ('pTM_ipTM', y_metric),       
                ('ipAE', y_metric),            
                ('pLDDT_peptide', y_metric),  
                ('pLDDT_peptide', 'pTM_ipTM'),
                ('pLDDT_peptide', 'ipTM')     
            ]
             
            # Limits & Ticks
            axis_limits = {
                'ipSAE': (-0.01, 0.85), 'ipSAE_min': (-0.01, 0.85),
                'ipAE': (6, 30),
                'pLDDT_peptide': (25, 82),
                'ipTM': (0.25, 0.99), 'pTM_ipTM': (0.25, 0.99)
            }
            axis_ticks = {
                'ipSAE': {'major': np.arange(0, 0.9, 0.2), 'minor': np.arange(0, 0.9, 0.1)},
                'ipSAE_min': {'major': np.arange(0, 0.9, 0.2), 'minor': np.arange(0, 0.9, 0.1)},
                'ipAE': {'major': np.arange(6, 31, 6), 'minor': np.arange(6, 31, 3)},
                'pLDDT_peptide': {'major': np.arange(40, 81, 20), 'minor': np.arange(30, 81, 10)},
                'ipTM': {'major': np.arange(0.3, 1.1, 0.2), 'minor': np.arange(0.3, 1.1, 0.1)},
                'pTM_ipTM': {'major': np.arange(0.3, 1.1, 0.2), 'minor': np.arange(0.3, 1.1, 0.1)}
            }
             
            fig, axes = plt.subplots(3, 2, figsize=fig_size)
            axes = axes.flatten()
            stats_data = []
             
            for i, (yy, xx) in enumerate(pairs_config):
                if yy not in df.columns or xx not in df.columns:
                    continue
                     
                ax = axes[i]
                valid = df.dropna(subset=[xx, yy])
                x_data = pd.to_numeric(valid[xx], errors='coerce')
                y_data = pd.to_numeric(valid[yy], errors='coerce')
                 
                # Calc Stats
                slope, intercept, r2, r = get_stats(x_data, y_data)
                stats_data.append({'Pair': f'{yy} vs {xx}', 'R_squared': r2, 'Pearson_r': r})
                 
                # Plot
                ax.scatter(x_data, y_data, alpha=0.6, edgecolors='w', s=dot_size, color='blue')
                 
                # Line
                if len(x_data) > 1:
                    lx = np.array([x_data.min(), x_data.max()])
                    ax.plot(lx, slope * lx + intercept, color='red', linewidth=2, linestyle='--')
                 
                # Text
                txt = f'$\mathbf{{r = {r:.3f}}}$'
                ax.text(0.05, 0.96, txt, transform=ax.transAxes, fontsize=TEXT_SIZE, fontweight='bold',
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))
                 
                # Limits
                if xx in axis_limits: ax.set_xlim(axis_limits[xx])
                if yy in axis_limits: ax.set_ylim(axis_limits[yy])
                 
                ax.set_xlabel(xx, fontsize=LABEL_SIZE, fontweight='bold')
                ax.set_ylabel(yy, fontsize=LABEL_SIZE, fontweight='bold')
                 
                apply_custom_style(ax, axis_ticks.get(xx), axis_ticks.get(yy))
                 
            plt.tight_layout()
            pname = os.path.join(output_dir, f'{base_name}_plots_all_{suffix}_16_18.png')
            plt.savefig(pname)
            plt.close()
             
            # Save Stats
            sdf = pd.DataFrame(stats_data)
            sname = os.path.join(output_dir, f'{base_name}_regression_stats_{suffix}.csv')
            sdf.to_csv(sname, index=False)
            print(f"Generated {suffix} plots: {pname}")

def _plot_scatter_ipsae_colored(complexes, source_file, start=1, end=None):
    """
    Generates Scatter Plots colored by ipSAE/ipSAE_min
    """
    print("Generating ipSAE-colored Scatter Plots...")
    
    # Imports for specific font handling if needed, but we rely on simple name
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.linewidth'] = 1.5
    LABEL_SIZE = 21
    TICK_SIZE = 20
    LEGEND_SIZE = 14
    
    for complex_name in complexes:
        work_home = os.getcwd()
        data_dir = os.path.join(work_home, "af2_init_guess", "data", complex_name)
         
        target_csv = None
        patterns = [
            os.path.join(data_dir, f"*_with_ipSAEmin_*.csv"),
            os.path.join(data_dir, f"{complex_name}_merged_scores_pTM-ipTM_with_oriSubs.csv")
        ]
        for p in patterns:
            matches = glob.glob(p)
            for m in matches:
                try:
                    d = pd.read_csv(m)
                    if 'ipSAE' in d.columns or 'ipSAE_min' in d.columns:
                        target_csv = m
                        break
                except: pass
            if target_csv: break
             
        if not target_csv: continue
         
        output_dir = os.path.dirname(target_csv)
        try:
            df = pd.read_csv(target_csv)
             
            metrics = {}
            if 'ipAE' in df.columns: metrics['y'] = 'ipAE'
            elif 'pae_interaction' in df.columns: metrics['y'] = 'pae_interaction'
            else: metrics['y'] = df.columns[7] if len(df.columns) > 7 else None
             
            if 'pTM_ipTM' in df.columns: metrics['x'] = 'pTM_ipTM'
            else: metrics['x'] = None 
             
            if not metrics['x'] or not metrics['y']: continue

            df[metrics['x']] = pd.to_numeric(df[metrics['x']], errors='coerce')
            df[metrics['y']] = pd.to_numeric(df[metrics['y']], errors='coerce')
             
            plot_configs = []
            if "ipSAE" in df.columns:
                plot_configs.append({
                    "col": "ipSAE", "vmin": 0.2, "vmax": 0.8, "label": r"$\mathbf{ipSAE}$", "suffix": "vs_ipSAE"
                 })
            if "ipSAE_min" in df.columns:
                plot_configs.append({
                    "col": "ipSAE_min", "vmin": 0.05, "vmax": 0.35, "label": r"$\mathbf{ipSAE}\_\mathbf{min}$", "suffix": "vs_ipSAE_min"
                 })
                 
            if not plot_configs: continue
             
            def setup_axis(ax):
                ax.set_ylabel(r'$\mathbf{ipAE}$ ($\AA$)', fontsize=LABEL_SIZE)
                ax.set_xlabel(r'$\mathbf{0.2pTM+0.8ipTM}$', fontsize=LABEL_SIZE)
                ax.tick_params(axis='both', which='major', labelsize=TICK_SIZE, length=6, width=1.7)
                ax.tick_params(axis='x', which='minor', length=4, width=1.5)
                for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontweight('bold')
                ax.yaxis.set_major_locator(ticker.MultipleLocator(6))
                ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
                ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
                ax.set_xlim(0.35, 0.96)
                ax.set_ylim(5.7, 30)
                for spine in ax.spines.values(): spine.set_linewidth(1.5)
                ax.set_title(f'{complex_name}', fontsize=10)

            size = 100
            size_star = 600
            star_edge = 4
            is_combined = "with_oriSubs" in target_csv or "combined" in target_csv

            for cfg in plot_configs:
                col = cfg['col']
                df[col] = pd.to_numeric(df[col], errors='coerce')
                valid_df = df.dropna(subset=[metrics['x'], metrics['y'], col])
                if valid_df.empty: continue
                 
                # Plot 2: Colored
                fig, ax = plt.subplots(figsize=(7, 6))
                cmap = plt.get_cmap('turbo_r')
                norm = plt.Normalize(vmin=cfg['vmin'], vmax=cfg['vmax'])
                 
                if is_combined:
                    ax.scatter(x=valid_df[metrics['x']][:-1], y=valid_df[metrics['y']][:-1], 
                                c=valid_df[col][:-1], cmap=cmap, norm=norm, s=size)
                     
                    last_idx = df.index[-1]
                    if last_idx in valid_df.index:
                        row = valid_df.loc[last_idx]
                        ax.scatter(x=row[metrics['x']], y=row[metrics['y']], 
                                    marker='*', s=size_star, c=[row[col]], cmap=cmap, norm=norm,
                                    edgecolor='black', linewidths=star_edge, zorder=10, label='Parental')
                        ax.legend(loc='lower left', fontsize=LEGEND_SIZE)
                else:
                    ax.scatter(x=valid_df[metrics['x']], y=valid_df[metrics['y']], 
                                c=valid_df[col], cmap=cmap, norm=norm, s=size)
                 
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                cbar = plt.colorbar(sm, ax=ax)
                cbar.set_label(cfg['label'], fontsize=LABEL_SIZE)
                cbar.ax.tick_params(labelsize=19)
                for l in cbar.ax.get_yticklabels(): l.set_fontweight('bold')
                 
                setup_axis(ax)
                figname = os.path.join(output_dir, f"{complex_name}_scatter_pTM-ipTM_vs_pae-inter_{cfg['suffix']}_withOriSubs.png")
                plt.savefig(figname, dpi=300, bbox_inches='tight')
                plt.close()
                 
                # Plot 3: Labelled
                fig, ax = plt.subplots(figsize=(7, 6))
                 
                if is_combined:
                    ax.scatter(x=valid_df[metrics['x']][:-1], y=valid_df[metrics['y']][:-1], 
                                c=valid_df[col][:-1], cmap=cmap, norm=norm, s=size)
                    last_idx = df.index[-1]
                    if last_idx in valid_df.index:
                        row = valid_df.loc[last_idx]
                        ax.scatter(x=row[metrics['x']], y=row[metrics['y']], 
                                    marker='*', s=size_star, c=[row[col]], cmap=cmap, norm=norm,
                                    edgecolor='black', linewidths=star_edge, zorder=10, label='Parental')
                        ax.legend(loc='lower left', fontsize=LEGEND_SIZE)
                else:
                    ax.scatter(x=valid_df[metrics['x']], y=valid_df[metrics['y']], 
                                c=valid_df[col], cmap=cmap, norm=norm, s=size)
                                
                cbar = plt.colorbar(sm, ax=ax)
                cbar.set_label(cfg['label'], fontsize=LABEL_SIZE)
                cbar.ax.tick_params(labelsize=19)
                for l in cbar.ax.get_yticklabels(): l.set_fontweight('bold')
                 
                setup_axis(ax)
                 
                lowest = valid_df.nsmallest(5, metrics['y'])
                for _, row in lowest.iterrows():
                    rid = str(row['id'])
                    label = rid
                    if "design" in rid:
                        idx = rid.find("design")
                        parts = rid[idx:].split("_")
                        if len(parts) > 1: label = f"d{parts[1]}"
                    plt.text(row[metrics['x']], row[metrics['y']], label, fontsize=10, fontweight='bold')
                     
                figname_lbl = os.path.join(output_dir, f"{complex_name}_scatter_pTM-ipTM_vs_pae-inter_{cfg['suffix']}_labelled_withOriSubs.png")
                plt.savefig(figname_lbl, dpi=300, bbox_inches='tight')
                plt.close()
                 
                print(f"Generated {cfg['suffix']} scatter plots for {complex_name}")
                 
        except Exception as e:
            print(f"Error generating ipSAE color plots for {complex_name}: {e}")
