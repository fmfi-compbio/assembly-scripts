FLYE_OPT = "-t 4"
BWA_OPT = "-t 4"
MINIMAP_OPT = "-t 4"
BLASTN_OPT = "-evalue 1e-5 -num_alignments 10000 -num_threads 4"
SAMPLE_SEED = 42
SCRIPT_PATH = "/opt/assembly-scripts"

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

# chromosome sizes
rule sizes:
    input:
        "{assembly}.fa"
    output:
        "{assembly}.sizes"
    shell:
        """
	faSize -detailed {input} > {output}
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

# sometimes snakemake does not correctly split the filename
# into two parts - help via special rule with ALN in the middle
rule asm_minimap2:
    input:
        fa1="{asm1}.fa", fa2="{asm2}.fa"
    output:
        "{asm1}-ALN-{asm2}.paf"
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

# output bed file for the first of the sequences to be compared by blast
rule blast2bed:
    input:
        "{name}.blast"
    output:
        "{name}.bed"
    shell:
        """
        perl -lane '$name="$F[2]:$F[4]-$F[5]:$F[10]"; print join("\t", @F[6,8,9], $name, @F[0,1])' {input} | sort -k1,1 -k2g > {output}
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
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
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

# nanopore aligned by minimap with paf format
rule Nanopore_minimap_paf:
    input:
         fa="{genome}.fa", fastq="{reads}-N.fastq.gz"
    output:
         paf="{genome}-{reads}-N.paf"
        
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# read connections between contigs
rule Nanopore_connections:
    input:
         fa="{genome}.fa", view="{genome}-{reads}-N.paf.view"
    output:
         "{genome}-{reads}-N.pairs"        
    shell:
        """
        sort -k 3,3 -k5,5g {input.view} > {output}.tmp1
        perl -lane 'die unless @F==11; if($o ne $F[2]) {{ $n=0; $o=$F[2]; }} if($F[1] eq "-") {{ @F[8,9]=@F[9,8]; }} print join("\t", $F[2]. "_" . $n, "R", @F[6,8,1]); print join("\t", $F[2]. "_" . ($n+1), "L", @F[6,9,1]); $n++;' {output}.tmp1 > {output}.tmp2
	sort -k1,1 {output}.tmp2 > {output}.tmp3
        join {output}.tmp3 {output}.tmp3  > {output}.tmp4
	perl -lane 'die unless @F==9; if($F[1] eq "L" && $F[5] eq "R") {{ print  join("\t", @F[0,2,3,4,6,7,8]); }}' {output}.tmp4 > {output}
	rm {output}.tmp[1234]
        """

# group and count read connections between contigs in windows
rule Nanopore_connections_size:
    input:
         "{name}-N.pairs"
    output:
         "{name}-N.pairs{size}kb"  
    shell:
       """
       perl -lane '$x={wildcards.size}*1000; $F[2]=int($F[2]/$x)*$x; $F[5]=int($F[5]/$x)*$x; print join("\t", @F[1..6])' {input} > {output}.tmp
       sort {output}.tmp | uniq -c | sort -k1gr > {output}
       rm {output}.tmp
       """

rule nanopore_sample:
    input:
         "{name}-N.fastq.gz"
    output:
         "{name}-sample-{minsize}kb-{frac}-N.fastq.gz"  
    shell:
       """
       zcat {input} | perl -ne 'BEGIN {{ srand({SAMPLE_SEED}); }} $s.=$_; if($.%4==0) {{ if(length($_)>{wildcards.minsize}*1000 && rand(1)<{wildcards.frac}) {{ print $s; }} $s=""; }}' | gzip -c > {output}
       """

# convert reads to fasta
rule fastq_gz_to_fasta:
    input:
         "{name}.fastq.gz"
    output:
         "{name}.fa"
    shell:
       """
       zcat {input} | perl -lane 'if($.%4==1) {{ $F[0]=~s/^@//; print ">$F[0]"; }} if($.%4==2) {{ print $_; }} ' > {output}
       """

# manual assembly
rule manual_assembly:
    input:
         fa="{name}-fasta.list", regions="{name}-regions.list"
    output:
         fa="{name}.fa", joins="{name}-joins.list"
    shell:
       """
       {SCRIPT_PATH}/manual_assembly.pl -j {output.joins} {input.fa} {input.regions} > {output.fa}
       """

rule self_aln:
    input:
        fa="{name}.fa"
    output:
        psl="{name}-self.psl", tab="{name}-self.tab"
    shell:
        """
        lastdb {output.psl}-tmp {input}
        lastal {output.psl}-tmp {input} -E1e-10 > {output.psl}.maf
        maf-convert psl {output.psl}.maf > {output.psl}.tmp
        perl -lane 'print unless $F[13] eq $F[9] && $F[0] >= 0.99*$F[10] && $F[0] >= 0.99*$F[14]' {output.psl}.tmp > {output.psl}
        maf-convert tab {output.psl}.maf > {output.tab}
        rm {output.psl}.maf {output.psl}.tmp
        rm {output.psl}-tmp.*
        """

rule chain:
    input:
        "{name}.psl"
    output:
        "{name}.chain"
    shell:
        """
	pslToChain {input} {output}
	"""

rule psl_png:
    input:
        "{name}.psl.tab"
    output:
        "{name}.psl.tab.png"
    shell:
        """
	last-dotplot {input} {output}
	"""

# e.g. _l100_id90 requires length 100 and id 90%, similarly 300_70
rule psl_filter:
    input:
        "{name}.psl"
    output:
        "{name}_l{length}_id{id}.psl"
    shell:
        """
	perl -lane 'die unless @F==21; next if $F[12]-$F[11]<{wildcards.length} || $F[16]-$F[15]<{wildcards.length} || ($F[0]+$F[2])<({wildcards.id}/100)*($F[0]+$F[1]+$F[2]+$F[3]+$F[5]+$F[7]); print' {input} > {output}
        """
