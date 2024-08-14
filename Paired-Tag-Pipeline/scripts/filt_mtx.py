import gzip
import os
import sys

def read_index(file_name):
    open_func = gzip.open if file_name.endswith('.gz') else open
    with open_func(file_name, 'rt') as f:
        return {i+1:line.strip() for i, line in enumerate(f)}


def read_matrix(filename, new_cell_id, old_cell_id, old_gene_id):
    gene_hits = {}
    new_cell_id_set = set(new_cell_id.values())
    with gzip.open(filename, 'rt') as f:
        next(f) # Skip first line
        next(f) # Skip second line
        next(f) # Skip third line
        for line in f:
            gene_idx, cell_idx, value = map(int, line.strip().split())
            cell_id = old_cell_id[cell_idx]
            gene_id = old_gene_id[gene_idx]
            if cell_id in new_cell_id_set:
                gene_hits.setdefault(gene_id, {})[cell_id] = value
    return gene_hits


            
def write_filtered(gene_hits, new_cell_id, old_gene_id, output):
    
    with gzip.open(f"{output}/barcodes.tsv.gz", 'wt') as out:
        for cell_idx in new_cell_id:
            out.write(f"{new_cell_id[cell_idx]}\n")
    
    new_gene_list = {}
    with gzip.open(f"{output}/features.tsv.gz", 'wt') as out:
        i = 0
        for gene_id in sorted(gene_hits):
            n_cells = len(gene_hits[gene_id])
            if n_cells >= 1:
                i += 1
                new_gene_list[gene_id] = i
                out.write(f"{gene_id}\n")

    n_genes = len(new_gene_list)
    n_cells = len(new_cell_id)
    n_values = sum(len(gene_hits[gene_id]) for gene_id in gene_hits)
    

    with gzip.open(f"{output}/matrix.mtx.gz", 'wt') as out:
        out.write("%%MatrixMarket matrix coordinate real general\n%\n")
        out.write(f"{n_genes} {n_cells} {n_values}\n")
      
        # this part can be optimized a lot
        for gene_id in new_gene_list:
            gene_idx = new_gene_list[gene_id]
            for cell_idx in new_cell_id:
                cell_id = new_cell_id[cell_idx]
                if cell_id in gene_hits[gene_id]:
                    value = gene_hits[gene_id][cell_id]
                    out.write(f"{gene_idx} {cell_idx} {value}\n")


if __name__ == "__main__":
    metadata = sys.argv[1] # "path-to-metadata"
    pre_mtx = sys.argv[2] # "path-to-raw-matrix"
    output = sys.argv[3] # "prefix-of-filtered-matrix"

    os.makedirs(output, exist_ok = True)

    new_cell_id = read_index(metadata)
    old_cell_id = read_index(f"{pre_mtx}/barcodes.tsv.gz")
    old_gene_id = read_index(f"{pre_mtx}/features.tsv.gz")
    gene_hits = read_matrix(f"{pre_mtx}/matrix.mtx.gz", new_cell_id, old_cell_id, old_gene_id)
    write_filtered(gene_hits, new_cell_id, old_gene_id, output)
