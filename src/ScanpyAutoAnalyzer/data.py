import pkg_resources
from os.path import isfile, join, exists
from pathlib import Path
import warnings
import shutil

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

def help():
    """Generates a formatted help string in Markdown table format for each .md file's help text."""
    # Locate the directory with .md files within the package
    lib_dir = pkg_resources.resource_filename('ScanpyAutoAnalyzer', 'data')
    lib_path = Path(lib_dir)

    # Get all .md files in the directory
    md_files = list(lib_path.glob("*.md"))

    # Header for the Markdown table
    help_md = "These are the available mardown file keys:\n"

    for md_file in md_files:
        # Use the file's stem as the title
        file_key = md_file.stem
        # Read the first non-code section as help text
        help_text = read_first_non_code_section(md_file)
        
        # Format each entry as a Markdown table row
        help_md += f"### {file_key}\n {help_text} \n"

    return help_md


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
