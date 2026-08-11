# importing default config
# overrinding - see ntoed in annot.yaml
import os
config_path = f"{workflow.basedir}/annot.yaml"
assert os.path.exists(config_path)
configfile: config_path

SCRIPT_PATH = workflow.basedir

# getting frequently used variables from config
GENETIC_CODE = config["GENETIC_CODE"]
MAX_INTRON = config["MAX_INTRON"]

rule rna_trim:
    input:
        fq1="{rnaseq}_1.fastq.gz",  fq2="{rnaseq}_2.fastq.gz"
    output:
        t1="{rnaseq}_trimmed_1.fastq.gz", t2="{rnaseq}_trimmed_2.fastq.gz",
        log="{rnaseq}_trimmed_1.fastq.log", json="{rnaseq}_trimmed_1.fastq.json",
        html="{rnaseq}_trimmed_1.fastq.html"
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} -o {output.t1} -O {output.t2} -j {output.json} -h {output.html} &>  {output.log}
        """

rule trinity_paired:
    input:
        fq1="{rnaseq}_trimmed_1.fastq.gz", fq2="{rnaseq}_trimmed_2.fastq.gz"
    output:
        fa="{rnaseq}_tr.fa", map="{rnaseq}_tr.gene_trans_map"
    params:
        tmp = "{rnaseq}_tmp_trinityp", log="{rnaseq}_tr_paired.log"
    shell:
         """
        Trinity --version > {params.log}
        Trinity {config[TRINITY_OPT]} {config[TRINITY_OPT_PAIRS]} --seqType fq --left {input.fq1} --right {input.fq2} --full_cleanup --output {params.tmp}
        mv {params.tmp}.Trinity.fasta {output.fa}
        mv {params.tmp}.Trinity.fasta.gene_trans_map {output.map}
        """

rule transcripts_blat:
    input:
       tr="{rnaseq}_tr.fa", genome="genome.fa"
    output:
       psl="{rnaseq}_tr.psl", log="{rnaseq}_tr.psl.log"
    params:
       tmp_psl="{rnaseq}_tr.tmp_psl"
    shell:
        """
        blat -maxIntron={MAX_INTRON} {input.genome} {input.tr} {params.tmp_psl}
        pslCDnaFilter {config[CDNAFILTER_OPT]} {params.tmp_psl} {output.psl} 2> {output.log}
        rm {params.tmp_psl}
        """

# intermediate step in transcript orfs
# files tr-exons.gtf and tr-CDS.gtf can be deleted
rule transcripts_exons:
    input:
       "{rnaseq}_tr.psl"
    output:
       "{rnaseq}_tr-exons.gtf"
    params:
       tmp_bed="{rnaseq}_tr.tmp.bed",
       tmp2_bed="{rnaseq}_tr.tmp2.bed",
       tmp_gp="{rnaseq}_tr.tmp.gp",
       tmp_gtf="{rnaseq}_tr.tmp.gtf",
    shell:
        """
	pslToBed {input} {params.tmp_bed}
        # make names unique
	perl -lane '$F[3] .= "_" . $.; print join("\\t", @F)' {params.tmp_bed} > {params.tmp2_bed}
        bedToGenePred {params.tmp2_bed} {params.tmp_gp}
        genePredToGtf file {params.tmp_gp} {params.tmp_gtf} 
        # we get strange gtf incl. start/stop...
        # keep only exons, remove strand, sort
        perl -F'"\\t"' -lane 'next unless $F[2] eq "exon"; $F[6]="."; print join("\\t", @F);' {params.tmp_gtf} | sort -k1,1 -k4,4g > {output}
        rm {params.tmp_bed} {params.tmp2_bed} {params.tmp_gp} {params.tmp_gtf}
        """

rule transcripts_orfs:
    input:
       "{rnaseq}_tr-CDS.gtf"
    output:
       "{rnaseq}_tr_orfs.gp"
    shell:
        """
	gtfToGenePred -genePredExt {input} {output}
        """

# no utrs, multiple orfs up to {value} of longest
# e.g. for value=0.75, 75% of longest
# identical ORFs clustered, all is renumbered
rule transcripts_orfs_multi:
    input:
       "{rnaseq}_tr-CDS_multi{value}.gtf"
    output:
       "{rnaseq}_tr_orfs_multi{value}.gp"
    shell:
        """
	# skip UTR
	perl -lane 'print if $F[2]=~/^(CDS|start_codon|stop_codon)$/' {input} > {output}.tmp.gtf
	# convert to gene pred
	gtfToGenePred -genePredExt {output}.tmp.gtf {output}.tmp.gp
	perl -lane '$F[0]="XYZ"; $F[11]="XYZ"; print join("\t", @F)' {output}.tmp.gp > {output}.tmp2.gp
	sort -k2,2 -k4,4g -k3 {output}.tmp2.gp | uniq > {output}.tmp3.gp
	perl -lane '$F[0]="orf$."; $F[11]="orf$."; print join("\t", @F)' {output}.tmp3.gp > {output}
	rm {output}.tmp.gp {output}.tmp2.gp {output}.tmp3.gp {output}.tmp.gtf
        """



rule bam_star:
    input:
       fq1="{rnaseq}_trimmed_1.fastq.gz", fq2="{rnaseq}_trimmed_2.fastq.gz",genome="genome.fa"
    output:
       bam="{rnaseq}_trimmed.bam"
    params:
       index="{rnaseq}_trimmed.bam-tmp",
       prefix="{rnaseq}_trimmed.bam-"
    shell:
        """
        mkdir {params.index}
        STAR {config[STAR_OPT]} --runMode genomeGenerate --genomeDir {params.index} --genomeFastaFiles {input.genome}  --genomeSAindexNbases 11
    	STAR {config[STAR_OPT]} {config[STAR_OPT2]} --genomeDir {params.index} --alignIntronMax {config[MAX_INTRON]} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.prefix} e --readFilesCommand zcat
    	rm -r {params.index}
    	rm {params.prefix}Log.progress.out
        # -F 260 uses only primary alignments
    	samtools view -b -F 260 -m {config[STAR_MIN_MATCH]} {params.prefix}Aligned.out.sam | samtools sort > {output.bam}
    	rm {params.prefix}Aligned.out.sam 
    	samtools index {output.bam}
        """

rule genome_sizes:
    input:
        "{genome}.fa"
    output:
        "{genome}.sizes"
    shell:
        """
        faSize -detailed {input} > {output}
        """

# create coverage from bam file (for rnaseq coverage)
rule cov_bedgraph:
    input:
        sizes="genome.sizes", bam="{aln}.bam"
    output:
        "{aln}.bedgraph"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -g {input.sizes} -bga -split > {output}
        """

# bigwig
rule cov_bg:
    input:
        sizes="genome.sizes", bedgraph="{aln}.bedgraph"
    output:
        "{aln}.bw"
    shell:
        """
        bedGraphToBigWig {input.bedgraph} {input.sizes} {output}
        """

rule repeats_fungi:
    input:
        "genome.fa"
    output:
        out="repeatMaskerFungi.out"
    params:
        tmp="tmp-rmf-genome.fa"
    shell:
        """
        cp {input} {params.tmp}
        chmod u+w {params.tmp}
        RepeatMasker -pa 4 -species fungi -xsmall {params.tmp}
        mv {params.tmp}.out {output}
        mv {params.tmp}.tbl {output}.tbl
        rm {params.tmp} {params.tmp}.cat.gz {params.tmp}.masked
        """

rule trna:
    input:
        "genome.fa"
    output:
        bed="trnascan.bed", out="trnascan.out"
    params:
        tmp="tm"
    shell:
        """
        tRNAscan-SE {input} > {output.out}
        {SCRIPT_PATH}/tRNAscan-SEtoBED.py < {output.out} > {output.bed}
        """

rule rfam:
    input:
        "genome.fa"
    output:
        bed="rfam.bed", out="rfam.out", tblout="rfam.tblout"
    shell:
        """
        cmscan {config[RFAM_OPT]} --rfam --cut_ga --nohmmonly --tblout {output.tblout} /extdata/Rfam/12.3/Rfam.cm {input} > {output.out}
        {SCRIPT_PATH}/cmscanToBed.py < {output.tblout} > {output.bed}
        """

rule hints:
    input:
        "transcripts.psl"
    output:
        "transcripts.hints.gff"
    shell:
        """
        sort -k14,14 -k16,16g {input} > {output}.tmp.psl
        {config[AUGUSTUS_DIR]}/scripts/blat2hints.pl --in={output}.tmp.psl --out={output}
        rm {output}.tmp.psl
        """

rule au_tr:
    input:
        fa="genome.fa", hints="transcripts.hints.gff", cfg="au-{cfg}.cfg"
    output:
        gtf="au-{cfg}-tr.orig.gtf"
    shell:
        """
        export SP=`head -n 1 {input.cfg}`
        export DIR=`head -n 2 {input.cfg} | tail -n 1`
        echo "SP:$SP  DIR:$DIR"
        augustus {config[AUGUSTUS_OPT]} --uniqueGeneId=true --AUGUSTUS_CONFIG_PATH=$DIR --species=$SP --hintsfile={input.hints} --extrinsicCfgFile={config[AUGUSTUS_DIR]}/config/extrinsic/extrinsic.ME.cfg {input.fa} > {output}
        """

rule au_gp:
    input:
        "au-{name}.orig.gtf"
    output:
        "au-{name,[a-zA-Z0-9_-]+}.gp"
    shell:
        """
        perl -lane 'print unless /^#/ || $F[1] ne "AUGUSTUS" || $F[2] eq "gene" || $F[2] eq "transcript"' {input} > {output}.tmp.gtf
        gtfToGenePred -genePredExt {output}.tmp.gtf {output}
        rm {output}.tmp.gtf
        """

rule au_gtf:
    input:
        gtf="au-{name}.gp"
    output:
        gtf="au-{name}.gtf"
    shell:
        """
        genePredToGtf -honorCdsStat file {input} {output}
        """

rule gtf2prot:
    input:
        gtf="{name}.gtf", fa="genome.fa"
    output:
        prot="{name}-prot.fa", cdna="{name}-cdna.fa"
    shell:
        """
        {SCRIPT_PATH}/gtf2transcript.pl -c -s -S -g {GENETIC_CODE} {input.fa} {input.gtf} -10 0 {output.cdna} {output.prot}
        """

# gtf should be sorted in the same order as genome.fa
rule exons2CDS:
    input:
        "{name}-exons.gtf"
    output:
        "{name}-CDS.gtf"
    shell:
        """
        {SCRIPT_PATH}/filter-gtf -p "INPUT {input} OUTPUT {output}" /opt/assembly-scripts/about/add_cds.about genome.fa
        """

# gtf should be sorted in the same order as genome.fa
# multiple orfs up to {value} of longest
# e.g. for value=0.75, 75% of longest
rule exons2CDS_multi:
    input:
        "{name}-exons.gtf"
    output:
        "{name}-CDS_multi{value}.gtf"
    shell:
        """
        {SCRIPT_PATH}/filter-gtf -p "INPUT {input} OUTPUT {output} FRACTION1 {wildcards.value} FRACTION2 0.25" {SCRIPT_PATH}/about/add_cds_multi.about genome.fa
        """


# list of supported transcripts for training
rule supported_by_transcripts:
    input:
        cdna="{name}-cdna.fa", tr="transcripts.fa"
    output:
        "{name}-supTr.list"
    shell:
        """
        blat -noHead -maxIntron=10 {input.tr} {input.cdna} {output}.tmp.psl
        # get only transcripts with 99% coverage by aln,
        #99% identity and no gaps
        #(gaps = potential bad introns or frameshifts)
        perl -lane '$m = $F[11]+($F[10]-$F[12]); $g=$F[0]+$F[2]; $gap=$F[5]+$F[7]; $b=$F[1]+$F[3]+$gap; print if $m<0.01*$F[10] && $b/($g+$b)<0.01 && $gap==0' {output}.tmp.psl  > {output}.tmp2.psl
        perl -lane 'print $F[9]' {output}.tmp2.psl | sort | uniq > {output}
        wc -l {output}
        rm {output}.tmp.psl {output}.tmp2.psl
        """

rule supported_by_other_prot:
    input:
        pred="{name}-prot.fa", other="other-prot.fa"
    output:
        "{name}-supProt.list"
    shell:
        """
        blat -noHead -prot {input.other} {input.pred} {output}.tmp.psl
	perl -lane 'print $F[9] if $F[0]>{config[AUGUSTUS_TRAIN_PROT_ID]}*$F[14] && $F[0]>{config[AUGUSTUS_TRAIN_PROT_ID]}*$F[10] && $F[8] eq "+"'  {output}.tmp.psl | sort -u > {output}
        wc -l {output}
        rm {output}.tmp.psl
        """

rule filter_supported:
    input:
        gtf="{name}.gtf", list="{name}-sup{type}.list"
    output:
        "{name}-sup{type,[a-zA-Z]+}.gtf"
    shell:
        """
	perl -lane 'print "transcript_id \\"$_\\";"' {input.list} > {output}.tmp.list
        grep -F -f {output}.tmp.list {input.gtf} > {output}
        rm {output}.tmp.list
        """

rule deoverlap_gtf:
    input:
        "{name}.gtf"
    output:
        "{name}-single.gtf"
    shell:
        """
        {SCRIPT_PATH}/deoverlap.pl {input} {output}
        """

# list of non-overlapping transcripts for training
rule list_for_training:
    input:
        "{name}-single.gtf"
    output:
        "{name}-single.list"
    shell:
        """
	perl -F'"\\t"' -lane 'if(/transcript_id "([^"]+)";/) {{ print $1; }}' {input} | sort -u > {output}
        """

# list of non-overlapping transcripts for training
rule augustus_training_gb:
    input:
        list="{name}-sup{type}-single.list",
        gtf="{name}.gtf",
        fa="genome.fa"
    output:
        "{name}-sup{type,[a-zA-Z]+}-single.train.gb"
    shell:
        """
        {config[AUGUSTUS_DIR]}/scripts/gtf2gff.pl  < {input.gtf} --out={output}.tmp.gff
        {config[AUGUSTUS_DIR]}/scripts/gff2gbSmallDNA.pl --good={input.list} {output}.tmp.gff genome.fa {config[AUGUSTUS_TRAIN_UTR]} {output}
        rm {output}.tmp.gff
"""

# augustus training
rule augustus_train:
    input:
        cfg="au-{cfg}.cfg"
    output:
        directory("au-{cfg}-train")
    shell:
        """
        export SP=`head -n 1 {input.cfg}`
        export DIR=`tail -n +2 {input.cfg} | head -n 1`
	export DATA=`tail -n +3 {input.cfg} | head -n 1`
        echo "SP:'$SP'  DIR:'$DIR'  DATA:'$DATA'"
	perl -le "die \\"empty line 3\\" unless length(\\"$DATA\\")>0"
	perl -le "die \\"$DATA does not exist\\" unless -r \\"$DATA\\""
	perl -le "die \\"bad dir $DIR\\" unless \\"$DIR\\" eq \\"{output}/config/\\""
	mkdir -p {output}/config
	cp -pr {config[AUGUSTUS_DIR]}/config/cgp/  {output}/config/
	cp -pr {config[AUGUSTUS_DIR]}/config/extrinsic/  {output}/config/
	cp -pr {config[AUGUSTUS_DIR]}/config/model/  {output}/config/
	cp -pr {config[AUGUSTUS_DIR]}/config/profile/  {output}/config/
	#cp -pr {config[AUGUSTUS_DIR]}/config/ {output}
        #rm -r {output}/species
        mkdir {output}/config/species
        cp -pr {config[AUGUSTUS_DIR]}/config/species/generic {output}/config/species/
        {config[AUGUSTUS_DIR]}/scripts/new_species.pl --AUGUSTUS_CONFIG_PATH=`pwd`/{output}/config --species=$SP
        perl -lne 's/stopCodonExcludedFromCDS false/stopCodonExcludedFromCDS true/; print' -i {output}/config/species/$SP/${{SP}}_parameters.cfg
        export AUGUSTUS_CONFIG_PATH=`pwd`/{output}/config ;{config[AUGUSTUS_DIR]}/bin/etraining --species=$SP $DATA &> {output}.log
	"""

# convert Uniprot orig fasta to version with just middle id
rule uniprot2simplify:
    input:
        fa="{name}UniprotOrig-prot.fa"
    output:
        fa="{name}Uniprot-prot.fa"
    shell:
        """
	perl -lne 'if(/>/) {{ s/>\\S+\\|(\\S+)\\|\\S+ />$1 / or die; }} print' {input.fa} > {output.fa}
	"""

# from uniprot fasta create a tsv file with id, gene_id, description
rule uniprot2tsv:
    input:
        fa="{name}Uniprot-prot.fa"
    output:
        tsv="{name}Uniprot-prot.tsv"
    shell:
        """
	perl -lne 'if(/>/) {{ /^>(\\S+) (.*) OS=.* GN=(\\S+)\\s.*$/ or die; print join("\\t", $1, $3, $2) }}' {input.fa} | sort -k1,1 > {output.tsv}
	"""

# convert uniprot fasta to fasta with gene names
#  commas in gene names will be changed to __
#  also non-unique gene names will get suffix __2, __3 etc 
rule uniprot2genes:
    input:
        fa="{name}Uniprot-prot.fa"
    output:
        fa="{name}Genes-prot.fa"
    shell:
        """
	perl -lne 'if(/>/) {{ s/>(.*) GN=(\\S+)\\s.*$/>$2/ or die; s/,/__/g; $c=2; $n=$_; while($seen{{$n}}) {{ $n=$_ . "__$c"; $c++; }}; $seen{{$n}}++; $_ = $n; }} print' {input.fa} > {output.fa}
	"""


# align proteins by miniprot
# gtf should be sorted in the same order as genome.fa
rule miniprot:
    input:
        fa="genome.fa", faa="{name}-prot.fa"
    output:
        gtf="{name}-prot.gtf"
    shell:
        """
	miniprot -G{MAX_INTRON} {config[MINIPROT_OPT]} --gtf {input.fa} {input.faa} > {output}
        """

# miniprot producing gff3 for later conversion to gp
rule miniprot2:
    input:
        fa="genome.fa", faa="{name}-prot.fa"
    output:
        gtf="{name}-prot.gff3"
    shell:
        """
	miniprot -G{MAX_INTRON} {config[MINIPROT_OPT]} --gff {input.fa} {input.faa} > {output}.tmp
	perl -lne 'next if /^##PAF/; if(/ID=(\\w+);/) {{ $o=$1; die "target $_" unless /Target=(\\S+)\\s/; $n=$1; if(exists $rev{{$n}}) {{ $i=2; while(exists $rev{{"${{n}}_c$i"}}) {{ $i++; }} $n="${{n}}_c$i"; }} $rev{{$n}}=$o; $fwd{{$o}} = $n; s/ID=$o/ID=$n/ or die "sub1 $_"; }} elsif (/Parent=(\\w+);/) {{ $o=$1; die "unknown $_" unless exists $fwd{{$o}}; $n=$fwd{{$o}}; s/Parent=$o/Parent=$n/ or die "sub2 $_"; }} s/Rank=/rank=/g; s/Identity=/identity=/g; s/Positive=/positive=/g; s/Frameshift=/frameshift=/g; s/StopCodon=/stopcodon=/g; s/Donor=/donor=/g; s/Acceptor=/acceptor=/g; print' {output}.tmp > {output}
	rm {output}.tmp
        """

rule miniprot_gp:
    input:
        "{name}-prot.gff3"
    output:
        "{name}-prot.gp"
    shell:
        """
	gff3ToGenePred {input} {output}
        """

# remove counters from multi-mapping proteins
# may not work if some protein names end in _c<number>
rule miniprot_browser_gp:
    input:
        "{name}-prot.gp"
    output:
        "{name}-prot-browser.gp"
    shell:
        """
	perl -F'"\\t"' -lane '$F[0]=~s/_c[0-9]+$//; print join("\\t", @F);' {input} > {output}
        """
     

# compare 2 protein fasta files by BLASTP
rule blastp:
    input:
        fa1="{seq1}.fa", fa2="{seq2}.fa"
    output:
        "{seq1}-BLASTP-{seq2}.blast"
    shell:
        """
        makeblastdb -in {input.fa1} -dbtype prot -out {output}.tmp
        blastp -db {output}.tmp -query {input.fa2} -task blastp -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" {config[BLASTP_OPT]} > {output}.tmp2
        rm {output}.tmp.p*
        perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  {output}.tmp2 | sort -k3,3 -k1gr > {output}
        rm {output}.tmp2
        """

# reduce blast results
rule reduce_blast:
    input: "{name}.blast"
    output: "{name}-reduced{id}-{best,[0-9.]+}.blast"
    shell:
      """
      # filter those with matches less than [id] fraction of length of both query and target
      perl -lane 'print if $F[0]>={wildcards.id}*$F[3] && $F[0]>={wildcards.id}*$F[7]' {input} > {output}.tmp
      # on both sides keep only prots within [best] fraction of best match 
      sort -k3,3 -k12gr {output}.tmp | perl -lane 'if($F[2] eq $o) {{ next unless $F[11]>={wildcards.best}*$m; }} else {{$o = $F[2]; $m=$F[11]; }} print ' > {output}.tmp2
      sort -k7,7 -k12gr {output}.tmp2 | perl -lane 'if($F[6] eq $o) {{ next unless $F[11]>={wildcards.best}*$m; }} else {{$o = $F[6]; $m=$F[11]; }} print ' > {output}
      wc -l {input} {output}.tmp {output}.tmp2 {output} > {output}.log
      rm {output}.tmp {output}.tmp2
      """

# reduce blast results
rule match1_blast:
    input: "{name1}-BLASTP-{name2}.tsv"
    output: "{name1}-BLASTP-{name2}-1match.tsv"
    shell:
      """
      perl -le 'print join("\\t", qw/query query_len target target_len matches num_aln_query num_aln_target/)' > {output}
      perl -lane '$n=scalar(@F); die unless $n>=4 && $n%2==0; $matches=($n-2)/2; @p=split ";", $F[3]; die unless @p==3; print join("\\t", @F[0,1,2], @p[0,1], $matches, $p[2])' {input} >> {output}
      """


# correspondence table
# for each target (column 3) count how many times in the table
# then create tsv sorted by query which contains all targets
# for query list its id and length, for each target its length, number of matches and how many times in the table
rule blast_table:
    input: "{name}.blast"
    output: "{name}.tsv"
    shell:
      """
      # get counts of targets
      perl -lane 'print $F[2]' {input} | sort | uniq -c | sort -k2 > {output}.tmp
      # sort by target to join counts in
      sort -k3,3 {input} > {output}.tmp2
      # join and select appropriate columns; keep score for now for sorting targets for each query
      join -1 2 -2 3 {output}.tmp {output}.tmp2 | perl -lane 'print join("\\t", @F[7,8,0], join(";", @F[4,2,1]), $F[12])' | sort -k1,1 -k5gr > {output}.tmp3
      perl -lane 'if($.==1 || $F[0] ne $o) {{ print "" if $.>1; printf "%s\\t\\%d", $F[0], $F[1]; }} printf "\\t%s\\t%s", $F[2], $F[3]; $o=$F[0]; END {{ print "" }} ' {output}.tmp3 > {output}
      rm {output}.tmp {output}.tmp2 {output}.tmp3
     """

# add gene names
rule match1_ext_blast:
    input: tsv="{name1}-BLASTP-{name2}Uniprot-prot{name3}1match.tsv", genes="{name2}Uniprot-prot.tsv"
    output: "{name1}-BLASTP-{name2}Uniprot-prot{name3}1match-ext.tsv"
    shell:
      """
      # sort our table for join with gene names, skip header
      tail -n +2 {input.tsv} | sort -k3,3 > {output}.tmp1
      # create header
      head -n 1 {input.tsv} | perl -lne 'print join("\\t", $_, "gene", "description")' > {output}
      # join our table with gene table, change order of columns and sort with query id
      join -1 3 -t $'\\t' {output}.tmp1 {input.genes} | perl -F'"\\t"' -lane 'die unless scalar @F==9; print join("\\t", @F[1,2,0,3..8])' | sort -k1,1 >> {output}
      rm {output}.tmp1
      """

# from paf.view alignment file in which special sequences of interest
# were used as db and genome as query
# create a bed file containing name of special sequence
# togther with % coverage and  %id; number of matches are used as score
rule paf_view_bed:
    input: "{name}.paf.view"
    output: "{name}.paf.bed"
    shell:
      """
      perl -lane '$cov=sprintf "%.1f", ($F[9]-$F[8])*100/$F[7]; $n=join("_", $F[6], "cov".$cov, "id".$F[10]); print join("\t", @F[2,4,5],$n,@F[0,1])' {input} | sort -k1,1 -k2g > {output}
      """


# find tandem repeats by tantan
# with a given max length
rule tantan:
    input: "genome.fa"
    output: bed="tantan-max{max}.bed", txt="tantan-max{max}.txt"
    shell:
      """
      tantan -f 4 -w {wildcards.max} {input}  > {output.txt}
      # the columns are 0:chr, 1:start, 2:end (BED-style),
      # 3:repeat_length, 4:number of repetitions, motif, occurrences
      perl -lane 'print join("\\t", @F[0..2], "rep$._$F[3]_$F[4]")' {output.txt} > {output.bed}
      """
