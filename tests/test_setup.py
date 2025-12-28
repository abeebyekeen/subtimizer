import pytest
import os
import shutil
from subtimizer.workflow.setup import setup_folders

def test_setup_initial_creation(tmp_path):
    """
    Test that setup_folders creates the expected directory structure.
    Uses pytest's tmp_path fixture to avoid messing up the actual filesystem.
    """
    # Create a dummy complexes file
    list_file = tmp_path / "complexes.dat"
    list_file.write_text("TestComplex_1\nTestComplex_2\n")
    
    # Change working directory to tmp_path so setup creates folders there
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        setup_folders(str(list_file), 'initial')
        
        # Verify
        assert os.path.isdir("TestComplex_1")
        assert os.path.isdir("TestComplex_1/AFcomplex")
        assert os.path.isdir("TestComplex_2")
    finally:
        os.chdir(cwd)

def test_setup_mpnn_creation(tmp_path):
    """Test mpnn folder creation."""
    list_file = tmp_path / "complexes.dat"
    list_file.write_text("TestComplex_1\n")
    
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        # Pre-create the initial structure
        os.makedirs("TestComplex_1/AFcomplex")
        
        setup_folders(str(list_file), 'mpnn')
        
        # Verify
        assert os.path.isdir("TestComplex_1/AFcomplex/mpnn_des")
    finally:
        os.chdir(cwd)
