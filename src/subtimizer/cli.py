import click
from subtimizer.workflow.setup import setup_folders
from subtimizer.workflow.folding import run_folding
from subtimizer.workflow.design import run_design
from subtimizer.workflow.analysis import analyze_recovery
from subtimizer.workflow.clustering import run_clustering
from subtimizer.workflow.preparation import prepare_for_folding
from subtimizer.workflow.pdb_utils import run_pdb_fix, fix_pdbs_in_dir
from subtimizer.workflow.validation import run_validation
from subtimizer.workflow.reporting import run_reporting

@click.group()
@click.version_option()
def main():
    """Subtimizer: Structure-Guided Design of Kinase Peptide Substrates.

    Yekeen et al. bioRxiv (2025)
    """
    pass

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file containing folder names (e.g., example_list_of_complexes.dat)')
@click.option('--type', '-t', type=click.Choice(['initial', 'mpnn', 'original']), default='initial', help='Type of setup to perform.')
def setup(file, type):
    """Set up folder structures for different stages of the workflow."""
    setup_folders(file, type)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum number of concurrent SLURM jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based) of complexes to process.')
@click.option('--end', type=int, default=None, help='End index (1-based) of complexes to process.')
@click.option('--mode', type=click.Choice(['batch', 'parallel']), default='batch', help='Execution mode: batch (default) or parallel.')
@click.option('--stage', type=click.Choice(['initial', 'validation']), default='initial', help='Folding stage: initial (Step 2) or validation (Step 7).')
def fold(file, max_jobs, start, end, mode, stage):
    """Run folding simulations (AlphaFold) for the complexes."""
    run_folding(file, max_jobs, start, end, mode, stage)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum number of concurrent SLURM jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based) of complexes to process.')
@click.option('--end', type=int, default=None, help='End index (1-based) of complexes to process.')
@click.option('--mode', type=click.Choice(['batch', 'parallel']), default='batch', help='Execution mode: batch (default) or parallel.')
def design(file, max_jobs, start, end, mode):
    """Run ProteinMPNN design."""
    run_design(file, max_jobs, start, end, mode)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
def analyze(file):
    """Analyze design results (Sequence Recovery)."""
    analyze_recovery(file)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum number of concurrent SLURM jobs.')
def cluster(file, max_jobs):
    """Cluster designed sequences using CD-HIT.
    
    Also attempts to summarize results to 'cluster_summary.dat'. 
    Note: Summarization might require jobs to finish first.
    """
    run_clustering(file, max_jobs)
    # Import locally to avoid circulars if any, but logic is in same file usually
    from subtimizer.workflow.clustering import summarize_clusters
    summarize_clusters(file)

if __name__ == '__main__':
    main()
@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
def prep_fold(file):
    """Prepare clustered sequences for validation folding."""
    prepare_for_folding(file)
@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum concurrent jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def fix_pdb(file, max_jobs, start, end):
    """Fix PDB files for initial guess validation."""
    run_pdb_fix(file, max_jobs, start, end)

@main.command(hidden=True)
@click.option('--dir', type=click.Path(exists=True), required=True)
def internal_fix_pdb(dir):
    """Internal command to run PDB fixing logic."""
    fix_pdbs_in_dir(dir)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum concurrent jobs.')
@click.option('--binder-path', envvar='DL_BINDER_DESIGN_PATH', help='Path to dl_binder_design/af2_initial_guess/predict.py')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def validate(file, max_jobs, binder_path, start, end):
    """Run AF2 Initial Guess."""
    run_validation(file, max_jobs, binder_path, start, end)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def report(file, start, end):
    """Generate final reports and swarm plots."""
    run_reporting(file, start, end)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--pae-cutoff', default="15", help='PAE cutoff (default: 15).')
@click.option('--dist-cutoff', default="15", help='Distance cutoff (default: 15).')
@click.option('--max-jobs', '-n', type=int, default=8, help='Maximum concurrent jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def ipsae(file, pae_cutoff, dist_cutoff, max_jobs, start, end):
    """Submit job for ipSAE evaluation."""
    from subtimizer.workflow.ipsae_runner import submit_ipsae_job
    submit_ipsae_job(file, pae_cutoff, dist_cutoff, max_jobs, start, end)

@main.command(hidden=True)
@click.option('--file', '-f', type=click.Path(exists=True), required=True)
@click.option('--pae-cutoff', default="15")
@click.option('--dist-cutoff', default="15")
@click.option('--max-jobs', '-n', type=int, default=8)
@click.option('--start', type=int, default=1)
@click.option('--end', type=int, default=None)
def internal_ipsae(file, pae_cutoff, dist_cutoff, max_jobs, start, end):
    """Internal command to run ipSAE logic."""
    from subtimizer.workflow.ipsae_runner import execute_ipsae_workflow
    execute_ipsae_workflow(file, pae_cutoff, dist_cutoff, max_jobs, start, end)
