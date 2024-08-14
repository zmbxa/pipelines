import os
import gzip
import sys

def merge_matrix(merge_list_file):
   
    with open(merge_list_file) as f:
        output = f.readline().strip()
        matrices = [line.strip().split() for line in f]

    os.makedirs(output, exist_ok = True)


    total = 0
    cell_list, feat_list, merge = {}, {}, {}
    for sample, prefix in matrices:
        print(f"Reading {sample}...")

        with gzip.open(f"{sample}/barcodes.tsv.gz", "rt") as f:
            cell_idx2id = {i: prefix + ":" + line.strip() for i, line in enumerate(f, start=1)}

        with gzip.open(f"{sample}/features.tsv.gz", "rt") as f:
            feat_idx2id = {i: line.strip() for i, line in enumerate(f, start=1)}

        with gzip.open(f"{sample}/matrix.mtx.gz", "rt") as f:
            next(f) # Skip first line
            next(f) # Skip second line
            next(f) # Skip third line
            for line in f:
                feat_idx, cell_idx, value = map(int, line.strip().split())
                cell_id = cell_idx2id[cell_idx]
                feat_id = feat_idx2id[feat_idx]
                merge.setdefault(feat_id, {})[cell_id] = value
                cell_list[cell_id] = 1
                feat_list[feat_id] = 1
                total += 1

    # Sort and write output files
    with gzip.open(f"{output}/feature.tsv.gz", 'wt') as f:
        for i, feat_id in enumerate(sorted(feat_list.keys()), start=1):
            feat_list[feat_id] = i
            f.write(f"{feat_id}\n")

    with gzip.open(f"{output}/barcodes.tsv.gz", 'wt') as f:
        for i, cell_id in enumerate(sorted(cell_list.keys()), start=1):
            cell_list[cell_id] = i
            f.write(f"{cell_id}\n")

    # Write matrix
    n_cell = len(cell_list)
    n_feat = len(feat_list)
    with gzip.open(f"{output}/matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n%\n")
        f.write(f"{n_feat} {n_cell} {total}\n")
        for feat_id in sorted(merge.keys()):
            for cell_id in sorted(merge[feat_id].keys()):
                f.write(f"{feat_list[feat_id]} {cell_list[cell_id]} {merge[feat_id][cell_id]}\n")



if __name__ == '__main__':
    merge_list_file = sys.argv[1]
    merge_matrix(merge_list_file)
