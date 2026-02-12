# Subtimizer

**A Computational Workflow for Structure-Guided Design of Potent and Selective Kinase Peptide Substrates**

[![DOI](https://zenodo.org/badge/doi/10.1021/acs.jcim.5c02430.svg?style=svg)](https://doi.org/10.1021/acs.jcim.5c02430)
[![PyPI version](https://badge.fury.io/py/subtimizer.svg)](https://badge.fury.io/py/subtimizer)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/abeebyekeen/subtimizer?style=flat-square)](https://github.com/abeebyekeen/subtimizer/releases)

---

## A. Contents

A. [Contents](#A-contents)   
B. [Overview](#B-overview)   
C. [Configuration](#C-configuration-important)   
D. [Prerequisites](#D-prerequisites)  
E. [Installation](#E-installation)  
F. [Usage](#F-usage)  
G. [Citation](#G-citation)

---

## B. Overview

**Subtimizer** provides an automated, structure-guided workflow for designing peptide substrates for kinases. It integrates **AlphaFold-Multimer** for structural modeling, **ProteinMPNN** for sequence design, and **AlphaFold2**-based interface evaluation of designed substrates.


## C. Configuration (Customizing SLURM Templates)

The workflow uses SLURM job scripts generated from templates. To customize these for your HPC environment (partition names, memory limits, modules):

1.  **Initialize local templates**:
    ```bash
    subtimizer init-templates
    ```
    This creates a `subtimizer_templates/` directory in your current folder with copies of all default scripts.

2.  **Edit the templates**:
    Open the files in `subtimizer_templates/` (e.g., `fold_template.sh`) and modify the `#SBATCH` directives or `module load` commands.

3.  **Run Subtimizer**:
    The tool will automatically detect and use your local templates instead of the package defaults.

## D. Prerequisites

1. [Install Anaconda or Miniconda](https://www.anaconda.com/products/distribution) or [Mamba](https://mamba.readthedocs.io)
2. [Install ColabFold](https://github.com/YoshitakaMo/localcolabfold)
> Add ColabFold to PATH using `export PATH="/PathTo/colabfold/localcolabfold/colabfold-conda/bin:$PATH"`
3. [Install ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
> Add ProteinMPNN to PATH using `export MPNN_PATH="/PathTo/ProteinMPNN/"`
4. [Get the code for af2_initial_guess](https://github.com/nrbennet/dl_binder_design)
> Add code to PATH using `export DL_BINDER_DESIGN_PATH="/PathTo/dl_binder_design/af2_initial_guess/predict.py"`
5. **SLURM**: This workflow is optimized for HPC environments using SLURM for job scheduling.

## E. Installation

### 1. Set Up Conda/Mamba Environments

#### Step A: Set Up Main Environment
Create a environment named `subtimizer_env` with Python>=3.9:

```bash
# Create the environment
mamba create -n subtimizer_env python=3.9 -y

# Activate the environment
mamba activate subtimizer_env
```

#### Step B: Set Up Worker Environments
AlphaFold and ProteinMPNN run in separate environments to avoid dependency conflicts. Create these environments using the provided YAML files in the repository root.

```bash
mamba env create -f af2_des_env.yaml
mamba env create -f mpnn_des_env.yaml
```

### 2. Install Subtimizer
While in the `subtimizer_env` environment, you can install the package via PyPI (recommended) or from source.

**Option A: Install from PyPI (Recommended)**
```bash
pip install subtimizer
```

**Option B: Install from Source**
Use this if you want to modify the code or templates.
```bash
git clone https://github.com/abeebyekeen/subtimizer.git
cd subtimizer
pip install .
```

### 3. Verify Installation

```bash
subtimizer --help
```


## F. Usage

The workflow is managed through the `subtimizer` command. Use `subtimizer --help` to see all available commands.

### Common Command Line Options

Most subtimizer commands (`fold`, `design`, `validate`, `fix-pdb`) accept the following options to control execution:

*   **`-n, --max-jobs <int>`**: Controls concurrency.
    *   Default is **4**. Increase this if you have more resources/GPUs available (e.g., `-n 8`).
    *   *Note: In `parallel` mode, this should match your SLURM script's layout.*
*   **`--start <int>` / `--end <int>`**: Process a subset of the list.
    *   Example: `--start 1 --end 10` (Processes items 1 through 10 in your input list).

### 1. Setup Project Structure

Change into your working directory:
```bash
cd examples
```

Initialize the directory structure for your kinase complexes.

**Input**: A file (e.g., `example_list_of_complexes.dat`) containing the list of folder names/complexes.
```text
AKT1_2akt1tide
ALK_axltide
SGK1_1akt1tide
TEC_srctide
```

For example, you can create this file using:
```bash
echo -e "AKT1_2akt1tide\nALK_axltide\nSGK1_1akt1tide\nTEC_srctide" > example_list_of_complexes.dat
```


Once you have your list of complexes, you can initialize the project structure using:
```bash
subtimizer setup --file example_list_of_complexes.dat --type initial
```
This creates the project directories and necessary subfolders that Subtimizer expects to find.

### 2. Run AlphaFold-Multimer
Launch AlphaFold-Multimer for the listed complexes.

>**Note**: This step expects a FASTA file (e.g., `AKT1_2akt1tide.fasta`) to exist inside each complex folder. See the `examples` folder for an example of how to prepare the FASTA files.

#### Option A: Batch Mode (Default)
Submits individual jobs for each complex using `fold_template.sh`.
```bash
subtimizer fold --file example_list_of_complexes.dat --max-jobs 4
```

#### Option B: Parallel Mode (Multi-GPU)
Submits a **single** job (`run_fold_parallel.sh`) that manages a pool of parallel tasks on a multi-GPU node.
```bash
subtimizer fold --file example_list_of_complexes.dat --mode parallel --max-jobs 4
```
*   `--max-jobs`: Number of parallel tasks (should match the number of GPUs requested in `fold_parallel_template.sh`).
*   `--start` / `--end`: Optionally specify explicit range of complexes to process on the list (e.g., `--start 1 --end 10`).

### 3. Run ProteinMPNN Design

Perform sequence design on the generated structures.

#### 3.1 Setup MPNN Design Folders and Configurations
```bash
subtimizer setup --file example_list_of_complexes.dat --type mpnn
```

>Edit the **`design_config.json`** file created in your working directory to customize `chains_to_design` (default: "B") or `fixed_positions` (default: "4") for specific complexes.

#### 3.2 Run ProteinMPNN Design

##### Option A: Batch Mode (Default)
Submits individual jobs using `design_template.sh`.
```bash
subtimizer design --file example_list_of_complexes.dat --max-jobs 4
```

##### Option B: Parallel Mode
Submits a single multi-sequence job using `design_parallel_template.sh`.
```bash
subtimizer design --file example_list_of_complexes.dat --mode parallel --max-jobs 4 --start 1 --end 10
```

### 4. Analyze Design Results

Analyze sequence recovery.

```bash
subtimizer analyze --file example_list_of_complexes.dat
```

> The `analyze` command generates:
> * Combined FASTA files (`all_design.fa`)
> * Sequence Logos (`*_seqlogo.png`)
> * Sequence Recovery Plots (`sequence_recovery_stripplot.png` and `.csv`)

### 5. Sequence Clustering

Cluster designed sequences to remove duplicates and generate a summary file `cluster_summary.dat`.

```bash
subtimizer cluster --file example_list_of_complexes.dat
```

### 6. Preparing kinase-peptide (designed) for folding

Prepare sequences for folding with AlphaFold-Multimer.

```bash
subtimizer prep-fold --file example_list_of_complexes.dat
```

### 7. Fold designed sequences with AF-Multimer

>Note: The version of proteinMPNN used in this work does not generate pdbs. Hence the need for post-design folding.

>However, with the [newer version/LigandMPNN](https://github.com/dauparas/LigandMPNN) which generates structures of designed sequences, this step may not be necessary.

##### Option A: Batch Mode (Default)
Run on a single node.
```bash
subtimizer fold --file example_list_of_complexes.dat --stage validation --max-jobs 4
```

##### Option B: Parallel Mode (Multi-GPU)
Distribute the folding of designed sequences across multiple GPUs on a single node.
```bash
subtimizer fold --file example_list_of_complexes.dat --stage validation --mode parallel --max-jobs 4
```

> **Tip for Multi-Node Parallelism**: To scale up to multiple nodes (e.g., 4 nodes), launch the parallel command 4 times with different ranges:
>1. `subtimizer fold ... --start 1 --end 2` (Node 1)
>2. `subtimizer fold ... --start 3 --end 4` (Node 2)
>... and so on.

>**Note**: This requires manually creating different SLURM jobs or running from different interactive sessions.

### 8. Prepare PDBs for AF2 initial guess

>**Note**: af2_init guess has two requirements for the input pdb:
>1. the binder (substrate) has to be the first chain
>2. no overlapping residue numbers between chains

```bash
subtimizer fix-pdb --file example_list_of_complexes.dat
```

### 9. Validation (AF2 Initial Guess)

Run AlphaFold with initial guess.

> **Configuration**: This step requires the path to the `af2_initial_guess` code (specifically `predict.py`). You can provide this path via the `--binder-path` argument or the `DL_BINDER_DESIGN_PATH` environment variable.

> You can set the environment variable as follows:
> ```bash
> export DL_BINDER_DESIGN_PATH="/path/to/dl_binder_design/af2_initial_guess/predict.py"
> ```

**Command**:
```bash
subtimizer validate --file example_list_of_complexes.dat --binder-path /path/to/dl_binder_design/af2_initial_guess/predict.py
```

### 10. Reporting

```bash
subtimizer report --file example_list_of_complexes.dat
```

> This `report` command generates final reports, including:
> * Merged score CSVs with weighted `pTM_ipTM` metric (`0.2*pTM + 0.8*ipTM`).
> * Swarm plots of validation metrics.
> * Data is copied to `af2_init_guess/data/` for easy access.

### 11. Workflow for Parental Substrates

To process the parental substrates, use the `setup --type original` command with the standard `subtimizer` workflow tools as described in the previous sections.

1.  **Setup**: Creates `original_subs` folder and prepares files.
    ```bash
    subtimizer setup --file example_list_of_complexes.dat --type original
    ```

2.  **Process**: Run commands on the parental substrates.
    ```bash
    cd original_subs
    # Fix PDBs
    subtimizer fix-pdb --file ../example_list_of_complexes.dat
    
    # Validation
    subtimizer validate --file ../example_list_of_complexes.dat --max-jobs 4
    
    # Reporting (Generates 'original' data)
    subtimizer report --file ../example_list_of_complexes.dat
    ```

3.  **Final Merge**: Return to the main directory and run/re-run report to combine results.
    ```bash
    cd ..
    # This automatically detects 'original_subs' and merges the data
    subtimizer report --file example_list_of_complexes.dat
    ```

### 12. ipSAE Evaluation

Perform a secondary interface evaluation on the folded structures using the ipSAE metric.

1. **Prerequisite**: Download [ipSAE](https://github.com/dunbracklab/IPSAE) and add it to your PATH:
   ```bash
   export PATH=$PATH:/path/to/ipSAE_directory
   ```

2. **Run ipSAE Calculation**:
   This command submits a SLURM job to calculate ipSAE metrics for all structures (designed and parental).
   ```bash
   subtimizer ipsae --file example_list_of_complexes.dat --max-jobs 8
   ```
   The `ipsae` command supports:
   *   list slicing: `subtimizer ipsae --file list.dat --start 1 --end 5`
   *   ipSAE cutoffs: `--pae-cutoff` (default: 15), `--dist-cutoff` (default: 15).

3. **Generate Final Reports**:
   Run the report command again to generate ipSAE-specific plots (Regression and Scatter plots).
   ```bash
   subtimizer report --file example_list_of_complexes.dat
   ```

---

## G. Citation

If you use Subtimizer in your work, please cite:

> **Yekeen A.A., Meyer C.J., McCoy M., Posner B., Westover K.D.** Subtimizer: Computational Workflow for Structure-Guided Design of Potent and Selective Kinase Peptide Substrates. *Journal of Chemical Information and Modeling* (2026). [https://doi.org/10.1021/acs.jcim.5c02430](https://doi.org/10.1021/acs.jcim.5c02430)




