#!/usr/bin/env python
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser(description="Convert SICER2 scoreisland files to BED format")
    parser.add_argument("input", nargs="+", help="Input scoreisland file(s)")
    parser.add_argument("-o", "--output-dir", help="Output directory (default: current directory)")
    args = parser.parse_args()

    output_dir = args.output_dir or "."

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for input_file in args.input:
        input_files = glob.glob(input_file)
        if not input_files:
            print(f"Error: Input file(s) '{input_file}' do(es) not exist")
            continue

        for file_path in input_files:
            if not os.path.isfile(file_path) or not os.access(file_path, os.R_OK):
                print(f"Error: Input file '{file_path}' does not exist or cannot be read")
                continue

            file_name = os.path.basename(file_path).split(".sco")[0]
            output_file = os.path.join(output_dir, f"{file_name}_peak.bed")

            with open(file_path, "r") as f_in, open(output_file, "w") as f_out:
                for line_num, line in enumerate(f_in, start=1):
                    fields = line.strip().split()
                    chrom, start, end, score = fields[0], fields[1], fields[2], fields[3]
                    peak_id = f"peaks{line_num}"
                    bed_line = f"{chrom}\t{start}\t{end}\t{peak_id}\t{score}\n"
                    f_out.write(bed_line)

            print(f"Conversion completed for file '{file_path}'. Output file: {output_file}")

if __name__ == "__main__":
    main()