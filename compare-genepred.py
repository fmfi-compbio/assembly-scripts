# Read two positional arguments by argparse: 
# input genepred file and chromosome sizes 
# and optional arguments with minimum overlap default 10 and expansion coefficient default 0.05
import argparse

parser = argparse.ArgumentParser(description='Compare genes within a genepred file')
parser.add_argument('genepred', help='Input genepred file')
parser.add_argument('chrom_sizes', help='Chromosome sizes file')
parser.add_argument('--min_overlap', type=int, default=10, help='Minimum overlap (default: 10)')
parser.add_argument('--expansion', type=float, default=0.05, help='Expansion coefficient (default: 0.05)')

args = parser.parse_args()

# Read the chromosome sizes file and store the sizes in a dictionary
chrom_sizes = {}
with open(args.chrom_sizes, 'r') as f:
    for line in f:
        chrom, size = line.strip().split('\t')
        chrom_sizes[chrom] = int(size)

# Read the genepred file and deduplicate genes based on the key (chrom, strand, start, end, exon count, exon coords)
# For each key keep one full record, set of prefixes, full prefixes and full ids
key2gene = dict()
with open(args.genepred, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        
        key = "\t".join(fields[1:3]+fields[5:10])
        extend_by = int((int(fields[6]) - int(fields[5])) * args.expansion)
        start2 = max(0, int(fields[5]) - extend_by)
        end2 = min(chrom_sizes[fields[1]], int(fields[6]) + extend_by)
        locus = f"{fields[1]}:{start2}-{end2}"

        gene = {
            'id': fields[0],
            'chrom': fields[1],
            'start': int(fields[5]),
            'end': int(fields[6]),
            'key': key,
            'locus': locus,
            'ids': {fields[0]},
            'overlaps': set()
        }
        if key in key2gene:
            key2gene[key]['ids'] |= {fields[0]}
        else:
            key2gene[key] = gene

# Find overlaps between genes, keep for each gene a list of overlapping ids (one per key) 
genes = list(key2gene.values())
genes.sort(key=lambda g: (g['chrom'], g['start'], g['end']))
current_chrom = None
current_genes = []
for gene in genes:
    if gene['chrom'] != current_chrom:
        current_chrom = gene['chrom']
        current_genes = []
    new_current_genes = []
    for other in current_genes:
        if other['end'] < gene['start']:
            continue
        overlap_start = max(gene['start'], other['start'])
        overlap_end = min(gene['end'], other['end'])
        overlap_length = overlap_end - overlap_start
        if overlap_length >= args.min_overlap:
            gene['overlaps'] |= {other['id']}
            other['overlaps'] |= {gene['id']}
            new_current_genes.append(other)
    current_genes = new_current_genes
    current_genes.append(gene)



def get_prefix_sets(ids, exclude1 = set(), exclude2 = set()):
    prefixes1 = set()
    prefixes2 = set()
    for id in ids:
        prefix_list = id.split(':')
        assert len(prefix_list) >= 3, "bad name format: " + id
        prefixes1.add(prefix_list[0])
        prefixes2.add(f"{prefix_list[0]}:{prefix_list[1]}")
    return (prefixes1.difference(exclude1), prefixes2.difference(exclude2))

def format_prefixes(prefixes):
    return " ".join(sorted(prefixes))

# deal with prefixes and print the results
for gene in genes:
    (pr1, pr2) = get_prefix_sets(gene['ids'])
    (overlap_pr1, overlap_pr2) = get_prefix_sets(gene['overlaps'], pr1, pr2)

    print("\t".join([gene['id'], 
                     str(len(pr1)), str(len(overlap_pr1)),
                     str(len(pr2)), str(len(overlap_pr2)),
                     str(len(gene['ids'])), str(len(gene['overlaps'])),
                     gene['locus'],
                     format_prefixes(pr1), format_prefixes(overlap_pr1),
                     format_prefixes(pr2), format_prefixes(overlap_pr2),
                     format_prefixes(gene['ids']), format_prefixes(gene['overlaps'])
                     ]))
