import pkg_resources
from os.path import isfile, join, exists
from pathlib import Path
import warnings
import shutil

import jupytext
import papermill as pm
import re
import ast

def run_with_parameters(name, output="./", parameters=None):
    """
    Converts one of the intern scripts to a Jupyter notebook, extracts parameters, checks for required values, 
    and runs the notebook using Papermill with the specified parameters.
    
    Parameters:
        name (str): one of the included info() elements.
        output (str): Path to save the executed notebook.
        parameters (dict): Dictionary of parameter values to pass to the notebook.
        
    Raises:
        ValueError: If any required parameters are missing from the `parameters` argument.
    """
    
    # Step 1: Convert markdown to notebook format

    md_file = getScript4( name )
    notebook = jupytext.read(md_file)

    # Ensure parameters dictionary is provided
    if parameters is None:
        raise ValueError("A dictionary of parameters is required but was not provided.")

    # Step 2: Find the `parameters` cell and extract variable names
    required_params = {}
    for cell in notebook.cells:
        if "parameters" in cell.metadata.get("tags", []):
            # Find all variable assignments in the parameters cell
            param_lines = cell.source.splitlines()
            for line in param_lines:
                # Updated regex to match arrays (lists) as well as simple values
                match = re.match(r"(\w+)\s*=\s*(\[[^\]]*\]|\S+)", line)
                if match:
                    var_name = match.group(1).strip()
                    default_value = match.group(2).strip()
                    # Safely evaluate the default value (using ast.literal_eval instead of eval)
                    try:
                        required_params[var_name] = ast.literal_eval(default_value)
                    except ValueError:
                        # If the default is a string or something unparseable, fallback to None
                        required_params[var_name] = None
    if len(required_params) == 0:
        raise ValueError( "The notebook does not provide any parameter section - lib error!" )
    # Step 3: Check for missing parameters and collect them
    missing_params = []
    for param, default_value in required_params.items():
        if param not in parameters:
            if default_value is not None:
                # Assign default value and warn
                parameters[param] = default_value
                warnings.warn(f"Parameter '{param}' is missing. Using default value: {default_value}.")
            else:
                missing_params.append(param)

    if missing_params:
        raise ValueError(
            f"Missing required parameters: {', '.join(missing_params)}. "
            "Please provide values for these parameters."
        )

    # Step 4: Save the notebook to a temporary file
    temp_notebook = join( output, "temp_notebook.ipynb" )
    jupytext.write(notebook, temp_notebook)
    output_notebook = join( output, f"{name}.ipynb" )
    # Step 5: Execute the notebook with Papermill and the provided parameters
    pm.execute_notebook(temp_notebook, output_notebook, parameters=parameters)

    print(f"Executed notebook saved as {output_notebook}")

# Usage example:
# parameters = {"param1": 20, "param2": "example_input"}
# run_markdown_with_parameters("example.md", parameters=parameters)

def copy_script_to_path(name, path):
    """Copies a script from the package resources to the specified path."""
    # Get the resource file path within the package
    libFile = getScript4(name)
    
    # Check if the file exists
    if not exists(libFile):
        warnings.warn(f"File {name} does not exist in the package.")
        return

    # Define the destination path where the file should be copied
    destination = join(path, f"{name}.md")
    
    # Check if the destination path exists
    destination_dir = Path(path)
    if not destination_dir.exists():
        warnings.warn(f"Destination path '{path}' does not exist. Creating it.")
        destination_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy the file
    shutil.copy(libFile, destination)
    print(f"File {name} copied to {destination}")

def getScript4(name="ExampleAnalysis"):
    """Gets the file path for a specified script within the package."""
    libFile = pkg_resources.resource_filename('ScanpyAutoAnalyzer', join('data', f"{name}.md"))
    if not exists(libFile):
        warnings.warn(f"Sorry - {name} is not one of the available options.")
    return libFile

def info_md():
    """Generates a formatted help string in Markdown table format for each .md file's help text."""
    # Locate the directory with .md files within the package
    lib_dir = pkg_resources.resource_filename('ScanpyAutoAnalyzer', 'data')
    lib_path = Path(lib_dir)

    # Get all .md files in the directory
    md_files = list(lib_path.glob("*.md"))

    # Header for the Markdown table
    info = "These are the available mardown file keys:\n"

    for md_file in md_files:
        # Use the file's stem as the title
        file_key = md_file.stem
        # Read the first non-code section as help text
        help_text = read_first_non_code_section(md_file)
        
        # Format each entry as a Markdown table row
        info += f"### {file_key}\n {help_text} \n"

    return info


def read_first_non_code_section(filepath):
    """Reads the first non-code, non-comment section of a markdown file."""
    first_section = []
    with_data = False  # To check if we've started gathering data
    in_code_block = False
    in_comment = False

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()

            # Detect start and end of code blocks
            if line.startswith("```"):
                if with_data:
                    break
                in_code_block = not in_code_block
                continue

            # Detect start and end of comment blocks
            if line.startswith("---"):
                if with_data:
                    break
                in_comment = not in_comment
                continue

            # Skip lines within code or comment blocks
            if in_code_block or in_comment:
                continue

            # Add non-code, non-comment lines to first_section
            if line:
                first_section.append(line)
                with_data = True

    # Join lines with line breaks to reconstruct the section
    return "\n".join(first_section)
