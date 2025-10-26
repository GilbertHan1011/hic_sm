#!/bin/bash

# Mustache loop detection script for Hi-C data
# Usage: mustache_loop_detection.sh -i <mcool_file> -r <resolution> -o <output_dir> -p <prefix> -t <threads> -pt <pvalue_threshold> -st <sparse_threshold> -m <mustache_path>

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            MCOOL_FILE="$2"
            shift 2
            ;;
        -r|--resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -pt|--pvalue_threshold)
            PVALUE_THRESHOLD="$2"
            shift 2
            ;;
        -st|--sparse_threshold)
            SPARSE_THRESHOLD="$2"
            shift 2
            ;;
        -m|--mustache_path)
            MUSTACHE_PATH="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 -i <mcool_file> -r <resolution> -o <output_dir> -p <prefix> -t <threads> -pt <pvalue_threshold> -st <sparse_threshold> -m <mustache_path>"
            echo "  -i, --input: Input mcool file"
            echo "  -r, --resolution: Resolution for loop detection"
            echo "  -o, --output: Output directory"
            echo "  -p, --prefix: Output file prefix"
            echo "  -t, --threads: Number of threads"
            echo "  -pt, --pvalue_threshold: P-value threshold for loop detection"
            echo "  -st, --sparse_threshold: Sparse threshold for loop detection"
            echo "  -m, --mustache_path: Path to mustache executable"
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$MCOOL_FILE" || -z "$RESOLUTION" || -z "$OUTPUT_DIR" || -z "$PREFIX" || -z "$THREADS" || -z "$PVALUE_THRESHOLD" || -z "$SPARSE_THRESHOLD" || -z "$MUSTACHE_PATH" ]]; then
    echo "Error: Missing required parameters"
    echo "Use -h or --help for usage information"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Set output file path
OUTPUT_FILE="${OUTPUT_DIR}/${PREFIX}_loops.tsv"

# Run mustache loop detection
echo "Running mustache loop detection..."
echo "Input: ${MCOOL_FILE}"
echo "Resolution: ${RESOLUTION}"
echo "Output: ${OUTPUT_FILE}"
echo "Threads: ${THREADS}"
echo "P-value threshold: ${PVALUE_THRESHOLD}"
echo "Sparse threshold: ${SPARSE_THRESHOLD}"

"${MUSTACHE_PATH}" \
    -f "${MCOOL_FILE}" \
    -r "${RESOLUTION}" \
    -p "${THREADS}" \
    -pt "${PVALUE_THRESHOLD}" \
    -st "${SPARSE_THRESHOLD}" \
    -o "${OUTPUT_FILE}"

# Check if the command was successful
if [[ $? -eq 0 ]]; then
    echo "Mustache loop detection completed successfully"
    echo "Output saved to: ${OUTPUT_FILE}"
else
    echo "Error: Mustache loop detection failed"
    exit 1
fi
