FLYE_OPT = "-t 4"
BWA_OPT = "-t 4"
BLASTN_OPT = "-evalue 1e-5 -num_alignments 10000 -num_threads 4"

# run flye assembler on nanopore data
rule flye:
    input:
        "{genome}-N.fastq.gz"
    output:
        fa="{genome}-flye.fa", log="{genome}-flye.fa.log"
    shell:
        """
        flye --version > {output.log}
        flye {FLYE_OPT} --nano-raw {input} --out-dir {output.fa}.tmp 2>> {output.log}
        mv {output.fa}.tmp/assembly.fasta {output.fa}
        mv {output.fa}.tmp/assembly_info.txt {output.fa}.assembly_info.txt
        mv {output.fa}.tmp/assembly_graph.gv {output.fa}.assembly_graph.gv
        mv {output.fa}.tmp/assembly_graph.gfa {output.fa}.assembly_graph.gfa
        gzip {output.fa}.assembly_graph.gfa
        """

# run alignments of illumina reads to an assembly
rule bwa:
    input:
        fa="{strain}-{assembly}.fa", il1="{strain}-I_1.fastq.gz", il2="{strain}-I_2.fastq.gz"
    output:
        bam="{strain}-{assembly}-I.bam", bai="{strain}-{assembly}-I.bam.bai",
    shell:
        """
      bwa index {input.fa} -p {output.bam}.tmp
      bwa mem {BWA_OPT} {output.bam}.tmp {input.il1} {input.il2} | samtools view -b - |  samtools sort - -o {output.bam}
      rm -f {output.bam}.tmp.amb {output.bam}.tmp.ann {output.bam}.tmp.bwt {output.bam}.tmp.pac {output.bam}.tmp.sa
      samtools index {output.bam}
       """

# pilon polishing with illumina bam file
rule pilon:
    input:
        fa="{assembly}.fa", bam="{assembly}-I.bam"
    output:
        fa="{assembly}p.fa", dir=directory("{assembly}-pilon")
    shell:
        """
        java -Xmx20G -jar /opt/broad/pilon-1.21.jar --genome {input.fa} --outdir {output.dir} --changes --tracks --frags {input.bam} > {output.fa}.log
        perl -lne 's/^(>.*)_pilon\s*$/$1/; print'  {output.dir}/pilon.fasta > {output.fa}
        """

# quast for any assembly
rule quast:
    input:
        "{assembly}.fa"
    output:
        "{assembly}.quast"
    shell:
        """
        quast.py -o {output}-tmp -t 1 {input}
        mv {output}-tmp/report.txt {output}
        """
        
# short summary of quast for logs
rule quastShort:
    input:
        "{assembly}.quast"
    output:
        "{assembly}.quastShort"
    shell:
        """
        perl -0777 -ne 'print "{wildcards.assembly} "; die unless /Total length \(>= 0 bp\)\s+(\d+)\s/; printf "Total %.1fM", $1/1e6; die unless /\# contigs \(>= 0 bp\)\s+(\d+)\s/; printf " contigs %d", $1; die unless /\# contigs \(>= 50000 bp\)\s+(\d+)\s/; printf " (%d >=50kb)", $1; die unless /Largest contig\s+(\d+)\s/; printf " longest %.0fkb", $1/1e3; die unless /N50\s+(\d+)\s/; printf " N50 %.0fkb", $1/1e3; die unless /\sGC \(\%\)\s+(\S+)\s/; printf " GC %.1f%%\n", $1; ' {input}  > {output}
        """

# compare two assemblies
rule asm_minimap:
    input:
        fa1="{asm1}.fa", fa2="{asm2}.fa"
    output:
        "{asm1}-{asm2}.paf"
    shell:
        """
        minimap2 -x asm20 -c {input.fa1} {input.fa2} > {output}
        """

# compare assembly to itself, skip trivial self-alignments
rule asm_minimap_self:
    input:
        fa="{asm}.fa"
    output:
        "{asm}-self.paf"
    shell:
        """
        minimap2 -x asm20 -c -p 0 {input.fa} {input.fa} | perl -lane 'print if $F[0] ne $F[5] || $F[4] ne "+" || $F[2] ne $F[7] || $F[3] ne $F[8]' > {output}
        """

# compare two fasta files by blastn, the first is treated as db
# the output is formatted similarly as paf.view files
rule blastn:
    input:
        fa1="{seq1}.fa", fa2="{seq2}.fa"
    output:
        "{seq1}-{seq2}.blast"
    shell:
        """
        makeblastdb -in {input.fa1} -dbtype nucl -out {output}.tmp
	blastn -db {output}.tmp -query {input.fa2} -task blastn -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" {BLASTN_OPT} > {output}.tmp2
        rm {output}.tmp.n*
	perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  {output}.tmp2 | sort -k3,3 -k1gr > {output}
	rm {output}.tmp2
        """


# create a pdf dotplot from a paf file
rule minimap_pdf:
    input:
        "{aln}.paf"
    output:
        "{aln}.paf.pdf"
    shell:
        """
        /usr/local/share/miniasm/miniasm/minidot -f 12 {input} | ps2pdf -dEPSCrop - {output}
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

# create a more readable version of a paf.gz file
rule minimap_gz_view:
    input:
        "{aln}.paf.gz"
    output:
        "{aln}.paf.view"
    shell:
        """
        zcat {input} | perl -lane '$id=sprintf("%.1f", $F[9]*100/$F[10]); print join("\t", @F[9,4,0..3,5..8],$id)' > {output}
        """

# nanopore aligned by minimap
rule Nanopore_minimap_bam:
    input:
         fa="{genome}.fa", fastq="{reads}-N.fastq.gz"
    output:
         bam="{genome}-{reads}-N.bam", bai="{genome}-{reads}-N.bam.bai"
        
    shell:
        """
        minimap2 -a -x map-ont --secondary=no -t 4 {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
    	samtools index {output.bam}
        """

# create coverage from bam file (for illumina / rnaseq coverage)
rule Nanopore_bedgraph:
    input:
         bam="{aln}.bam"
    output:
         cov="{aln}.bedgraph"
    shell:
        """
        bedgraph genomecov -ibam {input.bam} -bga -split {output.cov}
        """    