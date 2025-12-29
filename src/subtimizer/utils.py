import os
import shutil
import importlib.resources

TEMPLATE_DIR_NAME = "subtimizer_templates"

def get_template_content(template_filename):
    """
    Retrieves template content.
    Prioritizes local file in './subtimizer_templates/' if it exists.
    Falls back to the package default in 'subtimizer.templates'.
    """
    local_path = os.path.join(os.getcwd(), TEMPLATE_DIR_NAME, template_filename)
    
    if os.path.exists(local_path):
        print(f"Using custom template: {local_path}")
        with open(local_path, 'r') as f:
            return f.read()
    else:
        # Fallback to package
        try:
            return importlib.resources.files('subtimizer.templates').joinpath(template_filename).read_text()
        except FileNotFoundError:
            # Fallback for older python or if file strictly not found
            # try finding it via pkg_resources if importlib fails (shouldn't happen in py3.9+)
            raise FileNotFoundError(f"Template {template_filename} not found in package or local {TEMPLATE_DIR_NAME}.")

def copy_templates_to_local():
    """
    Copies all templates from the installed package to 'subtimizer_templates' in current dir.
    """
    dest_dir = os.path.join(os.getcwd(), TEMPLATE_DIR_NAME)
    os.makedirs(dest_dir, exist_ok=True)
    
    # List files in package
    # iterate over the resources object
    pkg_templates = importlib.resources.files('subtimizer.templates')
    
    count = 0
    for resource in pkg_templates.iterdir():
        if resource.is_file() and resource.name.endswith(".sh"):
            dest_path = os.path.join(dest_dir, resource.name)
            if not os.path.exists(dest_path):
                content = resource.read_text()
                with open(dest_path, 'w') as f:
                    f.write(content)
                print(f"Copied {resource.name} to {dest_dir}/")
                count += 1
            else:
                 print(f"Skipped {resource.name} (already exists)")
    
    if count > 0:
        print(f"\nSuccessfully copied {count} templates.")
        print("You can now edit these files to customize SLURM parameters.")
    else:
        print("\nNo new templates copied (files already existed).")
