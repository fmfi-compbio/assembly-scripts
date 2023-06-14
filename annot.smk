# HiSeq, MiSeq, palindrome mode for shorter inserts
TRIMMOMATIC_OPT = "ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10"
TRINITY_OPT = "--max_memory 40G --CPU 8"
BLAT_MAX_INTRON = "5000"
CDNAFILTER_OPT = "-minId=0.95 -minCover=0.75 -maxAligns=1"

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
        fa="{rnaseq}_tr_paired.fa", map="{rnaseq}_tr_paired.gene_trans_map",
        log="{rnaseq}_tr_paired.log"
    params:
        tmp = "{rnaseq}_tmp_trinityp"
    shell:
         """
        Trinity --version > {output.log}
        Trinity {TRINITY_OPT} --seqType fq --left {input.fq1} --right {input.fq2} --full_cleanup --output {params.tmp}
        mv {params.tmp}.Trinity.fasta {output.fa}
        mv {params.tmp}.Trinity.fasta.gene_trans_map {output.map}
        """

rule trinity_single:
    input:
        fq="{rnaseq}_trimmed_u.fastq.gz"
    output:
        fa="{rnaseq}_tr_single.fa", map="{rnaseq}_tr_single.gene_trans_map",
        log="{rnaseq}_tr_single.log"
    params:
        tmp = "{rnaseq}_tmp_trinitys"
    shell:
         """
        Trinity --version > {output.log}
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
        blat -maxIntron={BLAT_MAX_INTRON} {input.genome} {input.tr} {params.tmp_psl}
        pslCDnaFilter {CDNAFILTER_OPT} {params.tmp_psl} {output.psl} 2> {output.log}
        rm {params.tmp_psl}
        """