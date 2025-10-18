# Copied from https://github.com/mirnylab/distiller-sm/blob/master/_distiller_common.py
import os, pathlib
import numpy as np
import pandas as pd
import shlex

def argstring_to_dict(argstring):
    """
    Convert a command line argument string into a dictionary.
    
    Args:
        argstring (str): The command line argument string.
        
    Returns:
        dict: A dictionary with argument names as keys and their values.
        
    Note:
        - Arguments that start with '-' are considered keys.
        - If an argument has no value, it is set to True.
        - Values can be separated by spaces and will be joined into a single string.
        - If the final argument is a value without a preceding key, it will be included as a value for the last key (issue in case of specified input as last argument).
    """
    args = shlex.split(argstring)
    keys = np.where([arg.startswith('-') for arg in args])[0]
    vals = np.where([not arg.startswith('-') for arg in args])[0]
    args_arrs = [arr for arr in np.split(args, keys) if arr.size > 0]
    argdict = {str(arr[0]): (' '.join(arr[1:]) if arr.size > 1 else True) for arr in args_arrs }
    return argdict

def needs_downloading(fastq_files, side):
    if len(fastq_files) == 1 and fastq_files[0].startswith("sra:"):
        return True
    elif (
        fastq_files[side].startswith("sra:")
        or fastq_files[side].startswith("http://")
        or fastq_files[side].startswith("ftp://")
    ):
        return True
    else:
        return False

## Outputs a dictionary of the form {library: {run: [fastq1, fastq2]}}
def parse_fastq_folder(root_folder_path):
    library_run_fastqs = {}
    fastq_root_folder = pathlib.Path(root_folder_path)
    for fastq_file in sorted(fastq_root_folder.glob("**/*")):
        if not fastq_file.is_file():
            continue
        split_path = fastq_file.relative_to(fastq_root_folder).parts
        if len(split_path) != 3:
            raise Exception(
                "The fastq folder has a non-standard structure! "
                "Please, make sure that the fastq folders only contains a folder "
                "per library, each containing a folder per run, each _only_ containing "
                "two fastq(.gz) files."
            )
        library, run, fastq = split_path

        library_run_fastqs.setdefault(library, {}).setdefault(run, []).append(
            str(fastq_file.absolute())
        )

    return library_run_fastqs


def check_fastq_dict_structure(library_run_fastqs):
    for library, run_dict in library_run_fastqs.items():
        if not isinstance(run_dict, dict):
            return False
        for run, fastq_files in run_dict.items():
            if not isinstance(fastq_files, list) or (len(fastq_files) > 2):
                return False
    return True


def parse_fastq_dataframe(df):
    """
    Convert a pandas DataFrame to library_run_fastqs dictionary structure.
    
    Args:
        df (pd.DataFrame): DataFrame with columns: sample_id, lane, fastq1, fastq2
        
    Returns:
        dict: A nested dictionary with structure {library: {run: [fastq1, fastq2]}}
    """
    # Validate required columns
    required_columns = {'sample_id', 'lane', 'fastq1', 'fastq2'}
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        raise Exception(
            f"DataFrame is missing required columns: {missing}. "
            f"Required columns are: {required_columns}"
        )
    
    library_run_fastqs = {}
    
    # Iterate through DataFrame rows and build the nested dictionary
    for _, row in df.iterrows():
        library = str(row['sample_id'])
        run = str(row['lane'])
        fastq1 = str(row['fastq1'])
        fastq2 = str(row['fastq2'])
        
        # Initialize nested dictionaries if they don't exist
        if library not in library_run_fastqs:
            library_run_fastqs[library] = {}
        if run not in library_run_fastqs[library]:
            library_run_fastqs[library][run] = []
        
        # Add fastq files
        library_run_fastqs[library][run] = [fastq1, fastq2]
    
    return library_run_fastqs


def organize_fastqs(config):
    load_csv = config["input"]["load_csv"]
    if load_csv:
        fastq_csv = config["input"]["fastq_csv"]
        raw_reads_input = pd.read_csv(fastq_csv)
    else:
        raw_reads_input = config["input"]["raw_reads_paths"]
    
    # Case 1: String path to a folder
    if isinstance(raw_reads_input, str):
        library_run_fastqs = parse_fastq_folder(raw_reads_input)
    
    # Case 2: pandas DataFrame
    elif isinstance(raw_reads_input, pd.DataFrame):
        library_run_fastqs = parse_fastq_dataframe(raw_reads_input)
    
    # Case 3: Dictionary structure
    elif isinstance(raw_reads_input, dict):
        if not check_fastq_dict_structure(raw_reads_input):
            raise Exception(
                "An unknown format for library_fastqs! Please provide it as either "
                'a path to the folder structured as "library/run/fastqs", '
                "a dictionary specifying the project structure, or "
                "a pandas DataFrame with columns: sample_id, lane, fastq1, fastq2."
            )
        library_run_fastqs = raw_reads_input
    
    else:
        raise Exception(
            "An unknown format for library_fastqs! Please provide it as either "
            "a path to the folder with the structure library/run/fastqs, "
            "a dictionary specifying the project structure, or "
            "a pandas DataFrame with columns: sample_id, lane, fastq1, fastq2."
        )

    return library_run_fastqs
