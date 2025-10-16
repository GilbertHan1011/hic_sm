
from os import path
import re
from tempfile import TemporaryDirectory
from snakemake.shell import shell

# This script is called by Snakemake's 'script' directive.
# The 'snakemake' object is automatically available in this scope.

def basename_without_ext(file_path):
    """
    Calculates the base name of a file path, removing common FASTQ
    and compression extensions. This mimics FastQC's own naming logic.
    """
    base = path.basename(file_path)
    # A series of substitutions to robustly find the base name.
    base = re.sub(r"\.fastq\.gz$", "", base)
    base = re.sub(r"\.fq\.gz$", "", base)
    base = re.sub(r"\.fastq$", "", base)
    base = re.sub(r"\.fq$", "", base)
    base = re.sub(r"\.bam$", "", base)
    base = re.sub(r"\.sam$", "", base)
    base = re.sub(r"\.txt$", "", base)
    base = re.sub(r"\.gz$", "", base)
    base = re.sub(r"\.bz2$", "", base)
    return base

# Using a temporary directory is a best practice to prevent race conditions
# where multiple FastQC jobs might try to write to the same directory at once.
with TemporaryDirectory() as tempdir:
    # Execute the FastQC command using snakemake.shell
    shell(
        "fastqc"
        " --threads {snakemake.threads}"
        " --outdir {tempdir}"
        " {snakemake.params.extra}"
        " {snakemake.input[0]}"
        " &> {snakemake.log}"
    )

    # Determine the predictable output file names that FastQC creates.
    output_base = basename_without_ext(snakemake.input[0])
    temp_html_path = path.join(tempdir, output_base + "_fastqc.html")
    temp_zip_path = path.join(tempdir, output_base + "_fastqc.zip")

    # Move the generated files from the temporary directory to the
    # final destination specified in the 'output' directive.
    shell("mv {temp_html_path} {snakemake.output.html}")
    shell("mv {temp_zip_path} {snakemake.output.zip}")
