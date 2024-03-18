import os
import argparse
import glob

parser = argparse.ArgumentParser(description="This script is for extract bam file using given barcodes files.")
parser.add_argument("-b","--bam",nargs=1,help="REQUIRED, bam file that you need to be separated, like ./passorted_bam.bam",dest="bam_file")
parser.add_argument("-f",nargs="+",help="REQUIRED, barcode files",dest="barcode_files")
parser.add_argument("-o","--output_dir",help="output directory, ./split_bam if not specified",default="./split_bam",dest="output_dir")
parser.add_argument("-t","--thread",help="number od thread,default is 10",default=10)
args = parser.parse_args()

if not args.bam_file:
  print("BAM file is required!")
  parser.print_usage()
  exit(1)
  
if not args.barcode_files:
  print("Barcodes file is required!")
  parser.print_usage()
  exit(1)
  
if not os.path.exists(args.output_dir):
  os.makedirs(args.output_dir)
  
barcode_files = []
for pattern in args.barcode_files:
    matches = glob.glob(pattern)
    barcode_files.extend(matches)

if len(barcode_files) == 0:
    print("No barcode files found!")
    exit(1)
  
for file in args.barcode_files:
  type=os.path.splitext(os.path.basename(file))[0]
  print(type)
  tmp_sam="tmp.sam"
  
  os.system(f"sambamba view -t {args.thread} {args.bam_file} -H > {tmp_sam}")
  os.system(f"sambamba view -t {args.thread} {args.bam_file} | LC_ALL=C grep -F -f {file} >> {tmp_sam}")
  os.system(f"sambamba view -t {args.thread} -S {tmp_sam} -f bam -o {args.output_dir}/{type}.bam")
  
  if os.path.isfile(f"{args.output_dir}/{type}.bam"):
    print(f"Splited {args.bam_file} to {args.output_dir}/{type}.bam for: {file}")
    
if os.path.exists(tmp_sam):
  os.remove(tmp_sam)

  
  
  
  
