import argparse
import pandas as pd
from pathlib import Path

def parse_stats_file(stats_file):
    """
    Reads a stats file (tab-separated). If a line has more than two columns, only the first two are used.
    Returns a dictionary mapping key to value.
    """
    stats_dict = {}
    with open(stats_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    key, value = parts[0], parts[1]
                    stats_dict[key] = value
    return stats_dict

def main():
    parser = argparse.ArgumentParser(description="Concatenate library stats files into a summary table.")
    parser.add_argument("-i", "--input", required=True, nargs="+", help="Input stats files")
    parser.add_argument("-n", "--names", required=True, nargs="+", help="List of library names, one for each input file, in order")
    parser.add_argument("-o", "--output", required=True, help="Output summary TSV file")
    parser.add_argument("-l", "--log", default=None, help="Optional log file")

    args = parser.parse_args()

    if len(args.input) != len(args.names):
        raise ValueError("The number of --input files must match the number of --names provided.")

    # Dictionary to store all data
    all_stats = {}

    for stats_file, library_name in zip(args.input, args.names):
        stats_dict = parse_stats_file(stats_file)
        all_stats[library_name] = stats_dict

    # Convert to DataFrame
    df = pd.DataFrame(all_stats)

    # Save to TSV
    df.to_csv(args.output, sep='\t', index=True)

    if args.log:
        with open(args.log, 'w') as log_file:
            log_file.write(f"Combined {len(all_stats)} library stats files\n")
            log_file.write(f"Output written to {args.output}\n")

if __name__ == "__main__":
    main()