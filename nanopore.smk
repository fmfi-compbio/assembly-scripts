MINIMAP_OPT = "-t 4"

# run flye assembler on nanopore data
rule bamTofastq:
    input:
        "{name}-unaligned.bam"
    output:
        "{name}.fastq.gz"
    shell:
        """
	samtools fastq -T qs {input} | gzip -c > {output}
        """

rule fastq_gz_to_fasta:
    input:
         "{name}.fastq.gz"
    output:
         "{name}.fa"
    shell:
       """
       zcat {input} | perl -lane 'if($.%4==1) {{ $F[0]=~s/^@//; print ">$F[0]"; }} if($.%4==2) {{ print $_; }} ' > {output}
       """

rule minimap2:
    input:
        fa1="{asm1}.fa", fa2="{asm2}.fa"
    output:
        "{asm1}-ALN-{asm2}.paf"
    shell:
        """
        minimap2 -x asm20 -c {input.fa1} {input.fa2} > {output}
        """

# create a more readable version of a paf file
rule minimap_view:
    input:
        "{aln}.paf"
    output:
        "{aln}.paf.view"
    shell:
        """
        perl -lane '$id=sprintf("%.1f", $F[9]*100/$F[10]); print join("\t", @F[9,4,0..3,5..8],$id)' {input} > {output}
        """

# nanopore aligned by minimap with paf format
rule Nanopore_minimap_paf:
    input:
         fa="{genome}.fa", fastq="{reads}.fastq.gz"
    output:
         paf="{genome}-MAP-{reads}.paf"
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# nanopore aligned to ref.fa by minimap with paf format
rule Nanopore_minimap_paf_ref:
    input:
         fa="ref.fa", fastq="{reads}.fastq.gz"
    output:
         paf="{reads}.paf"
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# fasta stats
rule fasta_stats:
    input:
         fa="{name}.fa"
    output:
         stats="{name}.stats"
    shell:
        """
        faSize {input} |  perl -lane 'BEGIN {{ print "total reads mean max median"; }} if(/^\d/) {{ @a=@F; }} elsif(/^Total/) {{ print join(" ",$a[0], $a[11], $F[3], $F[10], $F[13]) }}' > {output}
        """

rule mapped_list:
     input: "{name}.paf.view"
     output: "{name}.paf.list"
     shell:
       """
       perl -lane 'print $F[2]' {input} | sort -u > {output}
       """

rule mapped_fa:
     input:
       list="{name}.paf.list",
       fa="{name}.fa"
     output: "{name}-mapped.fa"
     shell:
       """
       faSomeRecords {input.fa} {input.list} {output}
       """
