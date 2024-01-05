# HiSeq, MiSeq, palindrome mode for shorter inserts
TRIMMOMATIC_OPT = "ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10"
TRINITY_OPT = "--max_memory 40G --CPU 8"
TRINITY_OPT_PAIRS = ""  # default non-stranded; see https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
MAX_INTRON = "5000"
CDNAFILTER_OPT = "-minId=0.95 -minCover=0.75 -globalNearBest=0"
STAR_OPT = "--runThreadN 4"
STAR_OPT2 = "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2"
STAR_MIN_MATCH = 30
SCRIPT_PATH = "/opt/assembly-scripts"
RFAM_OPT = "--cpu 6"
AUGUSTUS_DIR = "/usr/local/share/augustus-3.2.3/augustus-3.2.3"
AUGUSTUS_OPT = "--alternatives-from-evidence=false --alternatives-from-sampling=false"
AUGUSTUS_TRAIN_UTR = 1000
AUGUSTUS_TRAIN_PROT_ID = 0.7
GENETIC_CODE = 1
MINIPROT_OPT = "-t4 -j1"  # j1 is non-mammalian splice sites
BLASTP_OPT = "-evalue 1e-5 -num_alignments 10000 -num_threads 4"

rule rna_trim:
    input:
        fq1="{rnaseq}_1.fastq.gz",  fq2="{rnaseq}_2.fastq.gz"
    output:
        t1="{rnaseq}_trimmed_1.fastq.gz", t2="{rnaseq}_trimmed_2.fastq.gz",
        u="{rnaseq}_trimmed_u.fastq.gz", u2="{rnaseq}_trimmed_u2.fastq.gz",
        log="{rnaseq}_trimmed_1.fastq.log"
    params:
        t1="{rnaseq}_trimmed_1.fastq", t2="{rnaseq}_trimmed_2.fastq",
        u="{rnaseq}_trimmed_u.fastq", u2="{rnaseq}_trimmed_u2.fastq",
        temp1="{rnaseq}_tmp_1.fastq", temp2="{rnaseq}_tmp_2.fastq"
    shell:
        """
        zcat {input.fq1} > {params.temp1}
        zcat {input.fq2} > {params.temp2}
        trimmomatic PE {params.temp1} {params.temp2} {params.t1} {params.u} {params.t2} {params.u2} {TRIMMOMATIC_OPT} 2> {output.log}
        rm {params.temp1} {params.temp2}
        gzip {params.t1} {params.u} {params.t2} {params.u2}
        """

rule trinity_paired:
    input:
        fq1="{rnaseq}_trimmed_1.fastq.gz", fq2="{rnaseq}_trimmed_2.fastq.gz"
    output:
        fa="{rnaseq}_tr_paired.fa", map="{rnaseq}_tr_paired.gene_trans_map"
    params:
        tmp = "{rnaseq}_tmp_trinityp", log="{rnaseq}_tr_paired.log"
    shell:
         """
        Trinity --version > {params.log}
        Trinity {TRINITY_OPT} {TRINITY_OPT_PAIRS} --seqType fq --left {input.fq1} --right {input.fq2} --full_cleanup --output {params.tmp}
        mv {params.tmp}.Trinity.fasta {output.fa}
        mv {params.tmp}.Trinity.fasta.gene_trans_map {output.map}
        """

rule trinity_single:
    input:
        fq="{rnaseq}_trimmed_u.fastq.gz"
    output:
        fa="{rnaseq}_tr_single.fa", map="{rnaseq}_tr_single.gene_trans_map"
    params:
        tmp = "{rnaseq}_tmp_trinitys", log="{rnaseq}_tr_single.log"
    shell:
        """
        Trinity --version > {params.log}
        Trinity {TRINITY_OPT} --seqType fq --single {input.fq} --full_cleanup --output {params.tmp}
        mv {params.tmp}.Trinity.fasta {output.fa}
        mv {params.tmp}.Trinity.fasta.gene_trans_map {output.map}
        """

rule transcripts:
    input:
        fa1="{rnaseq}_tr_paired.fa", fa2="{rnaseq}_tr_single.fa"
    output:
        fa="{rnaseq}_tr.fa"
    shell:
        """
        perl -ne 's/>TRINITY/>TRINITYP/; print' {input.fa1} > {output.fa}
        perl -ne 's/>TRINITY/>TRINITYS/; print' {input.fa2} >> {output.fa}
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
        pslCDnaFilter {CDNAFILTER_OPT} {params.tmp_psl} {output.psl} 2> {output.log}
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
	perl -lane '$F[3] .= "_" . $.; print join("\t", @F)' {params.tmp_bed} > {params.tmp2_bed}
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



rule bam_star:
    input:
       fq1="{rnaseq}_trimmed_1.fastq.gz", fq2="{rnaseq}_trimmed_2.fastq.gz",fqu="{rnaseq}_trimmed_u.fastq.gz", genome="genome.fa"
    output:
       bam="{rnaseq}_trimmed.bam"
    params:
       index="{rnaseq}_trimmed.bam-tmp",
       prefix="{rnaseq}_trimmed.bam-"
    shell:
        """
        mkdir {params.index}
        STAR {STAR_OPT} --runMode genomeGenerate --genomeDir {params.index} --genomeFastaFiles {input.genome}  --genomeSAindexNbases 11
    	STAR {STAR_OPT} {STAR_OPT2} --genomeDir {params.index} --alignIntronMax {MAX_INTRON} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.prefix} e --readFilesCommand zcat
    	STAR {STAR_OPT} {STAR_OPT2} --genomeDir {params.index} --alignIntronMax {MAX_INTRON} --readFilesIn {input.fqu} --outFileNamePrefix {params.prefix}u- e --readFilesCommand zcat
    	rm -r {params.index}
    	rm {params.prefix}Log.progress.out {params.prefix}u-Log.progress.out
        # -F 260 uses only primary alignments
    	samtools view -b -F 260 -m {STAR_MIN_MATCH} {params.prefix}Aligned.out.sam | samtools sort > {params.prefix}1.bam
    	samtools view -b -F 260 -m {STAR_MIN_MATCH} {params.prefix}u-Aligned.out.sam | samtools sort > {params.prefix}2.bam
    	rm {params.prefix}Aligned.out.sam {params.prefix}u-Aligned.out.sam
    	samtools merge {output.bam} {params.prefix}1.bam {params.prefix}2.bam
    	rm {params.prefix}1.bam {params.prefix}2.bam
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
        cmscan {RFAM_OPT} --rfam --cut_ga --nohmmonly --tblout {output.tblout} /extdata/Rfam/12.3/Rfam.cm {input} > {output.out}
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
        {AUGUSTUS_DIR}/scripts/blat2hints.pl --in={output}.tmp.psl --out={output}
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
        augustus {AUGUSTUS_OPT} --uniqueGeneId=true --AUGUSTUS_CONFIG_PATH=$DIR --species=$SP --hintsfile={input.hints} --extrinsicCfgFile={AUGUSTUS_DIR}/config/extrinsic/extrinsic.ME.cfg {input.fa} > {output}
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
        /opt/assembly-scripts/filter-gtf -p "INPUT {input} OUTPUT {output}" /opt/assembly-scripts/about/add_cds.about genome.fa
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
	perl -lane 'print $F[9] if $F[0]>{AUGUSTUS_TRAIN_PROT_ID}*$F[14] && $F[0]>{AUGUSTUS_TRAIN_PROT_ID}*$F[10] && $F[8] eq "+"'  {output}.tmp.psl | sort -u > {output}
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
        {AUGUSTUS_DIR}/scripts/gtf2gff.pl  < {input.gtf} --out={output}.tmp.gff
        {AUGUSTUS_DIR}/scripts/gff2gbSmallDNA.pl --good={input.list} {output}.tmp.gff genome.fa {AUGUSTUS_TRAIN_UTR} {output}
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
	cp -pr {AUGUSTUS_DIR}/config/cgp/  {output}/config/
	cp -pr {AUGUSTUS_DIR}/config/extrinsic/  {output}/config/
	cp -pr {AUGUSTUS_DIR}/config/model/  {output}/config/
	cp -pr {AUGUSTUS_DIR}/config/profile/  {output}/config/
	#cp -pr {AUGUSTUS_DIR}/config/ {output}
        #rm -r {output}/species
        mkdir {output}/config/species
        cp -pr {AUGUSTUS_DIR}/config/species/generic {output}/config/species/
        {AUGUSTUS_DIR}/scripts/new_species.pl --AUGUSTUS_CONFIG_PATH=`pwd`/{output}/config --species=$SP
        perl -lne 's/stopCodonExcludedFromCDS false/stopCodonExcludedFromCDS true/; print' -i {output}/config/species/$SP/${{SP}}_parameters.cfg
        export AUGUSTUS_CONFIG_PATH=`pwd`/{output}/config ;{AUGUSTUS_DIR}/bin/etraining --species=$SP $DATA &> {output}.log
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
	miniprot -G{MAX_INTRON} {MINIPROT_OPT} --gtf {input.fa} {input.faa} > {output}
        """

rule miniprot2:
    input:
        fa="genome.fa", faa="{name}-prot.fa"
    output:
        gtf="{name}-prot.gff3"
    shell:
        """
	miniprot -G{MAX_INTRON} {MINIPROT_OPT} --gff {input.fa} {input.faa} > {output}.tmp
	perl -lne 'next if /^##PAF/; if(/ID=(\w+);/) {{ $o=$1; die "target $_" unless /Target=(\S+)\s/; $n=$1; if(exists $rev{{$n}}) {{ $i=2; while(exists $rev{{"${{n}}_$i"}}) {{ $i++; }} $n="${{n}}_$i"; }} $rev{{$n}}=$o; $fwd{{$o}} = $n; s/ID=$o/ID=$n/ or die "sub1 $_"; }} elsif (/Parent=(\w+);/) {{ $o=$1; die "unknown $_" unless exists $fwd{{$o}}; $n=$fwd{{$o}}; s/Parent=$o/Parent=$n/ or die "sub2 $_"; }} s/Rank=/rank=/g; s/Identity=/identity=/g; s/Positive=/positive=/g; s/Frameshift=/frameshift=/g; s/StopCodon=/stopcodon=/g; s/Donor=/donor=/g; s/Acceptor=/acceptor=/g; print' {output}.tmp > {output}
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

# compare 2 protein fasta files by BLASTP
rule blastp:
    input:
        fa1="{seq1}.fa", fa2="{seq2}.fa"
    output:
        "{seq1}-BLASTP-{seq2}.blast"
    shell:
        """
        makeblastdb -in {input.fa1} -dbtype prot -out {output}.tmp
        blastp -db {output}.tmp -query {input.fa2} -task blastp -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" {BLASTP_OPT} > {output}.tmp2
        rm {output}.tmp.p*
        perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  {output}.tmp2 | sort -k3,3 -k1gr > {output}
        rm {output}.tmp2
        """

# reduce blast results
rule reduce_blast:
    input: "{name}.blast"
    output: "{name}-reduced{id}-{best}.blast"
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


# correspondence table
# for each target (column 3) count how many times in the table
# then create tsv sorted by query which contains all targets
# for query list its id and length, for each target its length, number of matches and how many times in the table
rule blast_table:
    input: "{name}.blast"
    output: "{name}.tsv"
    shell:
      """
      perl -lane 'print $F[2]' {input} | sort | uniq -c | sort -k2 > {output}.tmp
      sort -k3,3 {input} > {output}.tmp2
      join -1 2 -2 3 {output}.tmp {output}.tmp2 | perl -lane 'print join("\t", @F[7,8,0], join(";", @F[4,2,1]))' | sort > {output}.tmp3
      perl -lane 'if($.==1 || $F[0] ne $o) {{ print "" if $.>1; printf "%s\t\%d", $F[0], $F[1]; }} printf "\t%s\t%s", $F[2], $F[3]; $o=$F[0]; END {{ print "" }} ' {output}.tmp3 > {output}
      rm {output}.tmp {output}.tmp2 {output}.tmp3
     """ 
