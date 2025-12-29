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
    """Subtimizer: Structure-Guided Design of Kinase Peptide Substrates. Yekeen et al. bioRxiv (2025).
    """
    pass

@main.command()
def init_templates():
    """Copy default templates to local directory for customization."""
    from subtimizer.utils import copy_templates_to_local
    copy_templates_to_local()

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file containing folder names (e.g., example_list_of_complexes.dat)')
@click.option('--type', '-t', type=click.Choice(['initial', 'mpnn', 'original']), default='initial', help='Type of setup to perform.')
def setup(file, type):
    """Set up folder structures (Step 1)."""
    setup_folders(file, type)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum number of concurrent SLURM jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
@click.option('--mode', type=click.Choice(['batch', 'parallel']), default='batch', help='Execution mode.')
@click.option('--stage', type=click.Choice(['initial', 'validation']), default='initial', help='Stage: initial (Structure Gen) or validation (Designed Seq Folding).')
def fold(file, max_jobs, start, end, mode, stage):
    """Run AlphaFold-Multimer (Step 2 or 7)."""
    run_folding(file, max_jobs, start, end, mode, stage)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum number of concurrent SLURM jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
@click.option('--mode', type=click.Choice(['batch', 'parallel']), default='batch', help='Execution mode.')
def design(file, max_jobs, start, end, mode):
    """Run ProteinMPNN sequence design (Step 3)."""
    run_design(file, max_jobs, start, end, mode)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
def analyze(file):
    """Analyze design results (Step 4)."""
    analyze_recovery(file)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum concurrent jobs.')
def cluster(file, max_jobs):
    """Cluster designed sequences (Step 5)."""
    run_clustering(file, max_jobs)
    from subtimizer.workflow.clustering import summarize_clusters
    summarize_clusters(file)

if __name__ == '__main__':
    main()

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
def prep_fold(file):
    """Prepare sequences for folding (Step 6)."""
    prepare_for_folding(file)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum concurrent jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def fix_pdb(file, max_jobs, start, end):
    """Fix and prepare input PDBs (Step 8)."""
    run_pdb_fix(file, max_jobs, start, end)

@main.command(hidden=True)
@click.option('--dir', type=click.Path(exists=True), required=True)
def internal_fix_pdb(dir):
    """Internal: Fix PDB logic."""
    fix_pdbs_in_dir(dir)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--max-jobs', '-n', type=int, default=4, help='Maximum concurrent jobs.')
@click.option('--binder-path', envvar='DL_BINDER_DESIGN_PATH', help='Path to predict.py')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def validate(file, max_jobs, binder_path, start, end):
    """Run AF2 Initial Guess (Step 9)."""
    run_validation(file, max_jobs, binder_path, start, end)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def report(file, start, end):
    """Generate reports and plots (Step 10)."""
    run_reporting(file, start, end)

@main.command()
@click.option('--file', '-f', type=click.Path(exists=True), required=True, help='Path to the list file.')
@click.option('--pae-cutoff', default="15", help='PAE cutoff.')
@click.option('--dist-cutoff', default="15", help='Distance cutoff.')
@click.option('--max-jobs', '-n', type=int, default=8, help='Maximum concurrent jobs.')
@click.option('--start', type=int, default=1, help='Start index (1-based).')
@click.option('--end', type=int, default=None, help='End index (1-based).')
def ipsae(file, pae_cutoff, dist_cutoff, max_jobs, start, end):
    """Run ipSAE evaluation (Step 12)."""
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
