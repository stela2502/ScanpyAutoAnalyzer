import pkg_resources
from os.path import isfile, join, isdir, exists
from pathlib import Path
import warnings
import shutil

def copy_script_to_path ( name, path ):
	"""Copies a script from the package resources to the specified path."""
    # Get the resource file path within the package

    libFile = getScript4( name )

    # Check if the file exists
    if not exists(libFile):
        warnings.warn(f"File {name} does not exist in the package.")
        return

    # Define the destination path where the file should be copied
    destination = join(path, f"{name}.md")

    # Copy the file
    shutil.copy(libFile, destination)
    print(f"File {name} copied to {destination}")


def getScript4 ( name="ExampleAnalysis" ):

	libFile = pkg_resources.resource_filename('ScanpyAutoAnalyzer',join( 'data', f"{name}.md" ))
	if not exists( libFile ):
		warnings.warn( f"Sorry - {name} is not one of the options I can give you:\n{help()}" )

	return ( libFile )


def help():
    """Generates a formatted help string for each .md file's help text."""
    # Locate the directory with .md files
    lib_file = pkg_resources.resource_filename('ScanpyAutoAnalyzer', '')  # Adjust path if needed
    lib_path = Path(lib_file)

    # Get all .md files in the directory
    md_files = list(lib_path.glob("*.md"))

    # List to collect each help text entry
    help_entries = []

    for md_file in md_files:
        # Use the file's stem as the title
        file_key = md_file.stem
        # Read the first non-code section as help text
        help_text = read_first_non_code_section(md_file)
        
        # Format each entry with a title and the help text
        entry = f"### {file_key}\n{help_text}\n"
        help_entries.append(entry)

    # Join all entries with a separator
    return "\n".join(help_entries)




def read_first_non_code_section(filepath):
    first_section = []
    in_code_block = False

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()

            # Detect start and end of code blocks
            if line.startswith("```"):
                in_code_block = not in_code_block
                continue

            # Skip lines within code blocks
            if in_code_block:
                continue

            # Add non-code lines to first_section
            if line:
                first_section.append(line)
            # Stop at the next empty line, indicating end of the section
            elif first_section:
                break

    # Join lines with line breaks to reconstruct the section
    return "\n".join(first_section)