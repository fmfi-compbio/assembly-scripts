## Creating UCSC browser files and a manual annotation support

* Create several genepred files with `.gp` extension (various gene predictions, ORFs found in RNA-seq transcripts and va protein homology)
* For each file select two-level prefix: first prefix puts it into a group, e.g. transcipts and second prefix specifies source within group, e.g. particular RNA-seq
* Create file `rename-genes.yaml` which contains among others section listig all .gp files without .gp but with RENAME-prefix1-prefix2 added, e.g.

```yaml
for_manual:
  - funannotateCall-RENAME-FA-C
  - au-2-tr-RENAME-AU-2T
  - au-1-tr-supTr-RENAME-AU-1TST
  - canPsy_25C1_tr_orfs_multi0.5-RENAME-TR-25
  - canPsy_8C1_tr_orfs_multi0.5-RENAME-TR-8
  - debHan-mapped_orfs-RENAME-L-DH
```
* The first of these files is considered **base** of the annotation and by default its genes will be kept, with some later marked for removal or other genes marked for addition.

* Also optionally create `browser_url` variable which will be prefixed before postions in the form `chr:start-end`
* Then run `rename-genes for_manual/combined.tsv` 
* This creates folder for_manual. For each gp it adds 2 versions - one with prefixes added to transcript and gene names and one with also UTR removed. The one with just renaming can be uploaded to the browser. The one without UTRs is used for comparison - transcripts with identical CDS set are considered equivalent even if their UTRs differ.
* It also concatentates all `NOUTR.gp` files in `combined.gp`
* Table `for_manual/combined.tsv` contains for each group of equivalent transcripts from comined.gp
  * 0:selected transcript id (first in the order of precedence given by order in yaml)
  * 1:num prefixes 2:num overlapping prefixes
  * 3:num extended prefixes 4:num overlapping prefixes
  * 5:num transcript ids 6:num overlapping tr ids
  * 7:locus
  * 8:prefixes (space separated) 9:overlapping prefixes
  * 10:extended prefixes 11:overlapping extended prefixes
  * 12:transcript ids 13:overlapping tr ids
  * Overlap requires at least 10 bases
* Table `for_manual/combined.tsv` can be manually filtered by various one-liners to select genes of interest. The result should be two files:
  * `for_manual/sel_col1.tsv` with lines corresponding to base transcripts to be investigated.
  * `for_manual/sel_col2.tsv` with lines corresponding to other transcripts to be investigated.

```bash
# manual filtration
# from base annot (prefix FA) take unsupported genes
perl -F'"\t"' -lane 'print if $F[8]=~/\bFA\b/ && $F[1]==1;' for_manual/combined.tsv > for_manual/sel_col1.tsv; wc -l for_manual/sel_col1.tsv
# from other annot take those that have at least 2 prefixes and do not overlap FA
perl -F'"\t"' -lane 'print if $F[8]!~/\bFA\b/ && $F[9]!~/\bFA\b/ && $F[1]>1;' for_manual/combined.tsv > for_manual/sel_col2.tsv; wc -l for_manual/sel_col2.tsv
```
* The next step is to create a table `rename-genes for_manual/for_manual.tsv`

## Manual annotation

* Upload table `for_manual/for_manual.tsv` to a spreadsheet and edit it
* Each transcript should have next to it `P`(plus, for keeping or adding it) or `M`(minus for deleting or not adding it)
* To replace base transcript by another one, add the new transcript to column 3 and `P` or `M` to column 4
* More rows with the first 4 columns filled can be added
* Empty fields are marked by a dot
* Columns with transcript ID (1 and 3) can contain multiple values separated by space or comma
* Download the result as `manual.tsv` in the current folder

Optional: manual transcripts
* If a new transcript should be created for some gene, add it to  a new .gp file in `for_manual` folder. Make sure it already satisfies no UTR and use of a special prefix. Also gene and transcript ids should be distinct form each other
* Use names of these transcripts in `manual.tsv`
* Add the name of the file to `rename-genes.yaml` as follows (for file `for_manual/manual-tr.gp`):
```yaml
for_manual_added:
  - manual-tr
```

To create the manual annotation when all is done, run `/opt/assembly-scripts/rename-genes manual.gp` and also checked the log for useful stats and potential problems

## Renaming

```bash
# edit filenames in rename-genes.yaml with config
/opt/assembly-scripts/rename-genes manual-prot.fa
# get a list of bad genes (stop codons and mtDNA predictions)
/opt/assembly-scripts/rename-genes omit-genes1.list
# do renaming
/opt/assembly-scripts/rename-genes genes1.gff3
# creates genes1.gff3  genes1.gp  genes1.gtf  genes1.id_map  genes1-prot.fa genes1-cdna.fa
```

Example of corresponding part of `rename-genes.yaml`:
```bash
# prefix of files with nuclear genes (.gtf and -prof.fa)
nucl: manual

# prefix of files with mtDNA genes
mt: mtDNA/mtdna-coding

# file with a list of transcripts to omit (file can be created by rename-genes)
omit: omit-genes1.list

# name of sequence with mtDNA
mtDNA: mtDNA

# prefix of output filenames
name: genes1

# prefix of gene names
prefix: CANPSY_p

# how many digits to use
justify: 4

# genetic code nuclear
nucl_code: 12

# genetic code mitochondrial
mt_code: 4
```
