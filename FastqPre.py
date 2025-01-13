import os
import argparse
from Bio import SeqIO

def parse_sequence_length(file_path):
    """Parse sequence_length.txt to get length ranges for each file."""
    length_ranges = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            file_name, length_range = line.split()
            if '~' in length_range:
                min_len, max_len = map(int, length_range.split('~'))
            else:
                min_len, max_len = 275, 375  # Default values if not specified properly
            length_ranges[file_name] = (min_len, max_len)
    return length_ranges

def filter_sequences(input_dir, output_dir, length_file):
    """Filter sequences in fastq files based on length criteria."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parse length ranges from the sequence_length.txt
    length_ranges = parse_sequence_length(length_file)

    # Default length range
    default_range = (275, 375)

    # Process each fastq file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, file_name)

            # Get length range for this file or use default
            min_len, max_len = length_ranges.get(file_name, default_range)

            # Filter sequences
            with open(output_path, 'w') as output_handle:
                for record in SeqIO.parse(input_path, 'fastq'):
                    seq_len = len(record.seq)
                    if min_len <= seq_len <= max_len:
                        SeqIO.write(record, output_handle, 'fastq')

def main():
    # Use argparse to parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter FASTQ files based on sequence length.")
    parser.add_argument("-l", "--length_file", required=True, help="Path to the sequence_length.txt file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory where filtered files will be saved.")

    args = parser.parse_args()

    # Use the current directory as the input directory
    input_dir = os.getcwd()

    filter_sequences(input_dir, args.output_dir, args.length_file)

if __name__ == "__main__":
    main()
