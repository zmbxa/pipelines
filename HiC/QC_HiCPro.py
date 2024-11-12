import os
import re
import argparse

parser = argparse.ArgumentParser(description="This script is for extract mapping metrics from HiC-Pro.")
parser.add_argument("-s","--sample",nargs="+",help="REQUIRED, sample folder name under hic_results, ",dest="samples")
parser.add_argument("-p","--path",nargs=1,help="path contains hic_results/ folder",default="./",dest="path")
parser.add_argument("-o","--out",help="output file name",default="./all_sample_qc.txt",dest="outfile")

args = parser.parse_args()
if not args.samples:
  print("BAM file is required!")
  parser.print_usage()
  exit(1)

        

out = open(args.outfile,'w')
out.write("\t".join(["Sample","reported_Pairs(uniq_map)%","multi_map%","dup%","map%","cisLong/all%","cisShort/all%","allValidPairs","cis/trans"])+"\n")
for sample in args.samples:
    pairStat = open(os.path.join(args.path[0],"hic_results", "stats", sample, f"{sample}.mpairstat"))
    for line in pairStat:
        words = line.strip().split("\t")
        if words[0] == "Unmapped_pairs":
            mapped = 100-float(words[2])
        elif words[0] == "Unique_paired_alignments":
            uniq = words[2]
        elif words[0] == "Multiple_pairs_alignments":
            multi = words[2]
    pairStat.close()

    vpStat = open(os.path.join(args.path[0],"hic_results", "stats", sample, f"{sample}_allValidPairs.mergestat"))
    for line in vpStat:
        words = line.strip().split("\t")
        if words[0] == "valid_interaction":
            VP_all = int(words[1])
        elif words[0] == "trans_interaction":
            trans = int(words[1])
        elif words[0] == "cis_interaction":
            c_t = int(words[1])/trans
        elif words[0] == "valid_interaction_rmdup":
            dup = (int(VP_all) - int(words[1]))*100/int(VP_all)
        elif words[0] == "cis_longRange":
            cis_long = int(words[1])*100/VP_all
        elif words[0] == "cis_shortRange":
            cis_short = int(words[1])*100/VP_all
    vpStat.close()
    out.write("\t".join([sample,(uniq),(multi),"%.2f"%(dup),"%.2f"%(mapped),"%.2f"%(cis_long),"%.2f"%(cis_short),str(VP_all),"%.2f"%(c_t)])+"\n")
    
out.close()
