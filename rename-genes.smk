configfile: "rename-genes.yaml"

SCRIPT_PATH = "/opt/assembly-scripts"
MAKER_PATH = "/usr/local/share/maker/maker/bin"

NUCL = config["nucl"]
MT = config["mt"]
NAME = config["name"]
OMIT = config["omit"]

rule gather_bad:
  input: 
    gtfN=NUCL + ".gtf",
    protN=NUCL + "-prot.fa",
    protM=MT + "-prot.fa"  
  output:
    OMIT
  params:
    mtDNA = config.get("mtDNA", "mtDNA"),
  shell:
    """
    # find genes with stop codons
    cat {input.protN} {input.protM} | \
    grep -E '>|\*' | grep -B 1 '\*' | grep '>' | perl -lne 's/^>// or die; print' > {output}.tmp
    # find mt genes in nucl
    perl -lane 'if ($F[0] eq "{params.mtDNA}") {{ die unless /transcript_id \\\"([^\\"]+)\\"/; print $1; }} ' {input.gtfN} >> {output}.tmp
    sort -u {output}.tmp > {output}
    rm {output}.tmp
    """

rule rename_genes:
  input:
    gtfN=NUCL + ".gtf",
    protN=NUCL + "-prot.fa",
    gtfM=MT + ".gtf",
    protM=MT + "-prot.fa",
    fa="genome.fa",
    omit = OMIT
  output:
    map = NAME + ".id_map",
    gff3=NAME + ".gff3",
    gtf=NAME + ".gtf",
    gp=NAME + ".gp",
    prot=NAME + "-prot.fa",
    cdna=NAME + "-cdna.fa"
  params:
    prefix = config["prefix"],
    justify = config.get("justify", "4"),
    mtDNA = config.get("mtDNA", "mtDNA"),
    root = NAME
  shell:
    """
    echo {input} {output}
    # check that mt only mitochondrial
    perl -lane 'die unless $F[0] eq "{params.mtDNA}"' {input.gtfM}
    # join nuclear and mt files    
    cat {input.gtfN} {input.gtfM} > {params.root}.tmp.gtf
    cat {input.protN} {input.protM} > {params.root}.tmp-prot.fa
    # omit bad genes
    faSomeRecords -exclude {params.root}.tmp-prot.fa {input.omit}  {params.root}.tmp2-prot.fa
    perl -lane 'print "transcript_id \\"$_\\""' {input.omit} > {params.root}.tmp.omit
    grep -F -f  {params.root}.tmp.omit -v {params.root}.tmp.gtf > {params.root}.tmp2.gtf

    # create genepred
    gtfToGenePred -genePredExt {params.root}.tmp2.gtf {params.root}.tmp.gp

    # add stop codons to CDS
    {SCRIPT_PATH}/stop2CDS.pl {params.root}.tmp2.gtf > {params.root}.tmp3.gtf

    # convert to gff3
    {MAKER_PATH}/genemark_gtf2gff3 {params.root}.tmp3.gtf > {params.root}.tmp.gff3
    # sort order for renaming
    faSize -detailed {input.fa} | perl -lane 'print $F[0], "\t", $.' >  {params.root}.tmp.sort

    # rename with locus tag
    {MAKER_PATH}/maker_map_ids --prefix {params.prefix} --suffix '-' --justify {params.justify} {params.root}.tmp.gff3 --sort_order {params.root}.tmp.sort > {output.map}.tmp
    # add 0 at the end
    perl -lane '$F[1]=~s/({params.prefix}\\d+)/${{1}}0/ or die "bad line $_"; print join("\t", @F)' {output.map}.tmp > {output.map}

    {MAKER_PATH}/map_gff_ids {output.map} {params.root}.tmp.gff3
    {MAKER_PATH}/map_fasta_ids {output.map} {params.root}.tmp2-prot.fa
    {MAKER_PATH}/map_data_ids -col 1 {output.map} {params.root}.tmp.gp # transcript id
    {MAKER_PATH}/map_data_ids -col 12 {output.map} {params.root}.tmp.gp # gene_id
    mv -v {params.root}.tmp.gff3 {output.gff3}
    mv -v {params.root}.tmp2-prot.fa {output.prot}
    mv -v {params.root}.tmp.gp {output.gp}
    genePredToGtf -honorCdsStat file {output.gp} {output.gtf}

    # create -cds file from gtf (badprot has bad genetic code)
    {SCRIPT_PATH}/gtf2transcript.pl -c -s -S {input.fa} {output.gtf} -10 0 {output.cdna} {params.root}.tmp.badprot.fa

    rm {params.root}.tmp.sort {params.root}.tmp.gtf {params.root}.tmp2.gtf {params.root}.tmp3.gtf {output.map}.tmp  {params.root}.tmp-prot.fa {params.root}.tmp.omit {params.root}.tmp.badprot.fa
    """

rule manual_annot:
  input: tsv="{name}-manual.tsv", gpNew="{name}-manual-new.gp", gpOld="{name}-manual-old.gp"
  output: gp="{name}-manual.gp", log="{name}-manual.log", gtf="{name}-manual.gtf"
  shell:
    """
    # check that tsv has no repeated IDs in columns with genes to be removed, added
    perl -F'"\\t"' -lane 'die "bad space in $_" if $F[0]=~/\\s/ || $F[1]=~/\\s/; print $F[0]; print $F[1]' {input.tsv} | sort | uniq -c | perl -lane 'next if $F[1] eq "KEEP" || $F[1] eq "DELETE" || $F[1] eq "-"; die "repeated id $_" if $F[0]>1'
    # to remove: from the first column except those marked as KEEP
    # (add tab after gene name to use in grep)
    perl -F'"\\t"' -lane 'print $F[0], "\\t" unless $F[1] eq "KEEP" || $F[0] eq "-"' {input.tsv} > {output.gp}.tmp-remove.list
    # to add: everything in second column except KEEP, DELETE
    perl -F'"\\t"' -lane 'next if $F[1] eq "KEEP" || $F[1] eq "DELETE"; print $F[1],"\\t"' {input.tsv} > {output.gp}.tmp-add.list
    # get genes to be added
    grep -F -f {output.gp}.tmp-add.list {input.gpNew} > {output.gp}.tmp-add.gp || true
    # get genes not to be removed
    grep -v -F -f {output.gp}.tmp-remove.list {input.gpOld} > {output.gp}.tmp-remain.gp || true
    # join the two files
    cat {output.gp}.tmp-add.gp {output.gp}.tmp-remain.gp | sort -k 2,2 -k4,4g > {output.gp}
    
    # make gtf
    genePredToGtf -honorCdsStat file {output.gp} {output.gtf}

    # check that transcript IDs are unique
    perl -lane 'print $F[0]' {output.gp} | sort | uniq -c | perl -lane 'die "duplicated id $_" if $F[0]>1'
    # make sure that transcript id and gene id different 
    perl -lane 'die "identical gen id and transcrit id $_" if $F[0] eq $F[11]' {output.gp}

    # count genes
    echo "Stats" > {output.log}
    echo -n " Remove " >> {output.log}
    perl -F'"\\t"' -lane 'print if $F[1] eq "REMOVE"' {input.tsv} | wc -l >> {output.log}
    echo -n " Replace " >> {output.log}
    perl -F'"\\t"' -lane 'print if $F[1] ne "REMOVE" && $F[0] ne "-"' {input.tsv} | wc -l >> {output.log}
    echo -n " Add " >> {output.log}
    perl -F'"\\t"' -lane 'print if $F[1] ne "REMOVE" && $F[0] eq "-"' {input.tsv} | wc -l >> {output.log}

    echo "Check that the counts are the same for added genes" >> {output.log}
    wc -l {output.gp}.tmp-add.list {output.gp}.tmp-add.gp >> {output.log}

    echo -n "Incomplete genes among added: " >> {output.log}
    grep -c incmpl {output.gp}.tmp-add.gp >> {output.log} || true

    echo "Check counts of remove vs remain vs old total:"  >> {output.log}
    wc -l {output.gp}.tmp-remove.list {output.gp}.tmp-remain.gp {input.gpOld} >> {output.log}
    
    echo -n "New gene count " >> {output.log}
    wc -l {output.gp} >> {output.log}

    echo "Overlaps longer than 30bp:" >> {output.log}

    overlapSelect -statsOutput -excludeSelf -overlapBases=30 {output.gp} {output.gp} stdout >> {output.log}

    rm {output.gp}.tmp-add.list {output.gp}.tmp-add.gp {output.gp}.tmp-remove.list {output.gp}.tmp-remain.gp
    head -n 30 {output.log}
    """