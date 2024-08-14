import csv
import os
import re

file_path='./fastq/'
filenames = [f for f in os.listdir(file_path) if f.endswith('.fastq.gz')]
sample_number = sorted(list(set([filename.split('_')[0] for filename in filenames])))

gap = eval(re.findall("\d+",sample_number[0])[0])+len(sample_number)//2-1


samples = []

## for each sample
for sample in sample_number:
    sample_name = sample
    sample_type = 'RNA' if int(re.findall("\d+",sample_name)[0]) > gap else 'DNA'
    r1_file = sample_type+"/fastq/"+sample+"_R1.fastq.gz"
    r2_file = sample_type+"/fastq/"+sample+"_R2.fastq.gz"

    # add sample information
    samples.append([sample_name, sample_type, r1_file, r2_file])

with open('sample.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['sampleID', 'Type', 'R1', 'R2'])
    writer.writerows(samples)

