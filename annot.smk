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
GENETIC_CODE = 1
MINIPROT_OPT = "-t4 -j1"  # j1 is non-mammalian splice sites

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
        "au-{name}.gp"
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
