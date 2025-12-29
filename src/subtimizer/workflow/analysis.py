import os
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import numpy as np

import subprocess

def analyze_recovery(file_path: str):
    """
    aggregates validation metrics (sequence recovery) and plots the results.
    """
    print(f"Analyzing results for complexes in {file_path}")
    
    with open(file_path, 'r') as f:
        complexes = [line.strip() for line in f if line.strip()]

    _combine_fastas(complexes)

    _generate_weblogo(complexes)

    _extract_recovery_data(complexes)

    _plot_recovery_strip(complexes)

def _combine_fastas(complexes):
    """
    Combines individual design FASTA files into one all_design.fa
    """
    print("Combining FASTA files...")
    for folder in complexes:
        seqs_dir = os.path.join(folder, "AFcomplex", "mpnn_out", "seqs")
        output_file = os.path.join(seqs_dir, "all_design.fa")
        
        if not os.path.exists(seqs_dir):
            print(f"Warning: {seqs_dir} not found.")
            continue
            
        fastas = glob.glob(os.path.join(seqs_dir, "*.fa"))
        if not fastas:
            print(f"No .fa files found in {seqs_dir}")
            continue
            
        with open(output_file, 'w') as outfile:
            for fasta in fastas:
                with open(fasta, 'r') as infile:
                    lines = infile.readlines()
                    if len(lines) > 2:
                        outfile.writelines(lines[2:])
                    else:
                        outfile.writelines(lines)
        print(f"Created {output_file}")

def _generate_weblogo(complexes):
    """
    Generate sequence logo using WebLogo
    """
    print("Generating WebLogos...")
    for folder in complexes:
        seqs_dir = os.path.join(folder, "AFcomplex", "mpnn_out", "seqs")
        input_file = os.path.join(seqs_dir, "all_design.fa")
        output_file = os.path.join(seqs_dir, f"{folder}_seqlogo.png")
        
        if not os.path.exists(input_file):
            continue
            
        cmd = [
            "weblogo",
            "-f", input_file,
            "-D", "fasta",
            "-o", output_file,
            "-F", "png_print",
            "-A", "protein",
            "-s", "large",
            "--errorbars", "NO",
            "--color-scheme", "chemistry",
            "--logo-font", "Arial-BoldMT",
            "-U", "probability"
        ]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"Generated WebLogo: {output_file}")
        except FileNotFoundError:
            print("Error: 'weblogo' command not found. Please install weblogo (pip install weblogo).")
        except subprocess.CalledProcessError as e:
            print(f"Error running weblogo for {folder}: {e}")

def _extract_recovery_data(complexes):
    """
    extract sequence recovery data from MPNN output
    """
    print("Extracting sequence recovery data...")
    for folder in complexes:
        design_file = os.path.join(folder, "AFcomplex", "mpnn_out", "seqs", "all_design.fa")
        output_file = os.path.join(folder, "AFcomplex", "mpnn_out", "seqs", "sec_recovery.dat")
        
        if not os.path.exists(design_file):
            print(f"Warning: Design file not found for {folder}")
            continue
            
        try:
            with open(output_file, "w") as dataout:
                with open(design_file, "r") as datain:
                    for line in datain:
                        if "=" in line:
                            parts = line.strip().split(",")
                            for part in parts:
                                if "seq_recovery" in part:
                                    val = part.split("=")[-1]
                                    dataout.write(f"{val}\n")
        except Exception as e:
            print(f"Error processing {folder}: {e}")

def _plot_recovery_strip(complexes):
    """
    Generates strip plots for sequence recovery.
    """
    print("Generating strip plots...")
    matplotlib.use('AGG')
    
    # Collect data
    secrec_data = {}
    for folder in complexes:
        path = os.path.join(folder, "AFcomplex", "mpnn_out", "seqs", "sec_recovery.dat")
        if os.path.exists(path):
            try:
                with open(path) as f:
                    secrec_data[folder] = [float(line.strip()) for line in f]
            except ValueError:
                print(f"Error reading floats from {path}")
    
    if not secrec_data:
        print("No data found to plot.")
        return

    try:
        data = pd.DataFrame.from_dict(secrec_data, orient="index").T
        long_data = data.melt(var_name="Kinase-peptide complex",
                              value_name="Fraction of peptide sequence recovered")
        
        plt.figure(figsize=(10, 6))
        
        # use colorcet palette if available, else default
        try:
            palette = sns.color_palette(cc.glasbey, n_colors=len(complexes))
        except:
            palette = "viridis"

        sns.stripplot(x="Kinase-peptide complex",
                      y="Fraction of peptide sequence recovered",
                      data=long_data, 
                      hue="Kinase-peptide complex",
                      alpha=0.5, 
                      palette=palette,
                      legend=False)

        plt.xticks(rotation=75)
        plt.yticks(np.arange(0, 1, 0.1))
        plt.ylim(-0.05, 1.05) 
        plt.tight_layout()
        
        output_png = "sequence_recovery_stripplot.png"
        plt.savefig(output_png, dpi=300)
        print(f"Plot saved to {output_png}")
        plt.clf()
        
        data.to_csv("sequence_recovery_data.csv")
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Error plotting: {e}")

