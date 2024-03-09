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

# run miniasm assembler on nanopore data
rule miniasm:
    input:
        "{genome}-N.fastq.gz"
    output:
        fa="{genome}-miniasm.fa"
    shell:
        """
        echo "minimap2 and miniasm versions:" > {output}.log
        minimap2 --version >> {output}.log
        miniasm -V >> {output}.log
        # Overlap reads
        minimap2 -x ava-ont {MINIMAP_OPT} {input} {input} | gzip -1 > {output}.paf.gz 
        # Layout
        miniasm -f {input} {output.fa}.paf.gz > {output.fa}.gfa 2>>{output}.log
        rm {output}.paf.gz
        perl -lane 'print ">$F[1]\n$F[2]" if $F[0] eq "S"' {output}.gfa > {output}
        gzip {output}.gfa
        """


        
# run alignments of illumina reads to an assembly
rule bwa:
    input:
        fa="{assembly}.fa", il1="{reads}-I_1.fastq.gz", il2="{reads}-I_2.fastq.gz"
    output:
        bam="{assembly}-MAP-{reads}-I.bam", bai="{assembly}-MAP-{reads}-I.bam.bai",
    shell:
        """
        bwa index {input.fa} -p {output.bam}.tmp
        bwa mem {BWA_OPT} {output.bam}.tmp {input.il1} {input.il2} | samtools view -b - |  samtools sort - -o {output.bam}
        rm -f {output.bam}.tmp.amb {output.bam}.tmp.ann {output.bam}.tmp.bwt {output.bam}.tmp.pac {output.bam}.tmp.sa
        samtools index {output.bam}
        """

# run alignments of illumina reads to an assembly for reads called illumina_[12].fastq.gz
rule bwa2:
    input:
        fa="{assembly}.fa", il1="illumina_1.fastq.gz", il2="illumina_2.fastq.gz"
    output:
        bam="{assembly}-I.bam", bai="{assembly}-I.bam.bai",
    shell:
        """
       bwa index {input.fa} -p {output.bam}.tmp
       bwa mem {BWA_OPT} {output.bam}.tmp {input.il1} {input.il2} | samtools view -b - |  samtools sort - -o {output.bam}
       rm -f {output.bam}.tmp.amb {output.bam}.tmp.ann {output.bam}.tmp.bwt {output.bam}.tmp.pac {output.bam}.tmp.sa
       samtools index {output.bam}
       """

# run alignments of illumina reads to an assembly for reads called illumina-sample_[12].fastq.gz
rule bwa_sample:
    input:
        fa="{assembly}.fa", il1="illumina-sample_1.fastq.gz", il2="illumina-sample_2.fastq.gz"
    output:
        bam="{assembly}-IS.bam", bai="{assembly}-IS.bam.bai",
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

rule pilon_unconfirmed:
     input:
        fa="{assembly}.fa", dir="{assembly}-pilon", sizes="{assembly}.sizes"
     output:
        "{assembly}-unconfirmed.bw"
     shell:
         """
         wigToBigWig {input.dir}/pilonUnconfirmed.wig {input.sizes} {output}
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

# compare two assemblies via minimap, note ALN in the middle
rule asm_minimap2:
    input:
        fa1="{asm1}.fa", fa2="{asm2}.fa"
    output:
        "{asm1}-ALN-{asm2}.paf"
    shell:
        """
        minimap2 -x asm20 -c {input.fa1} {input.fa2} > {output}
        """

# compare two assemblies via minimap, specify -p param after ALN
rule asm_minimap2_p:
    input:
        fa1="{asm1}.fa", fa2="{asm2}.fa"
    output:
        "{asm1}-ALN{p}P-{asm2}.paf"
    shell:
        """
        minimap2 -x asm20 -c -p {wildcards.p} {input.fa1} {input.fa2} > {output}
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
        "{seq1}-ALN-{seq2}.blast"
    shell:
        """
        makeblastdb -in {input.fa1} -dbtype nucl -out {output}.tmp
	blastn -db {output}.tmp -query {input.fa2} -task blastn -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" {BLASTN_OPT} > {output}.tmp2
        rm {output}.tmp.n*
	perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  {output}.tmp2 | sort -k3,3 -k1gr > {output}
	rm {output}.tmp2
        """

# output bed file for the first of the sequences compared by blast
rule blast2bed:
    input:
        "{name}.blast"
    output:
        "{name}.blast.bed"
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

# compare protein fasta file (db) and dna file by blastx
# the output is formatted similarly as paf.view files
rule blastx:
    input:
        fa1="{seq1}.fa", fa2="{seq2}.fa"
    output:
        "{seq1}-ALNX-{seq2}.blast"
    shell:
        """
        makeblastdb -in {input.fa1} -dbtype prot -out {output}.tmp
	blastx -db {output}.tmp -query {input.fa2} -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" {BLASTN_OPT} > {output}.tmp2
        rm {output}.tmp.p*
	perl -lane '$F[6]--; ($s,$e)=@F[2,3]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1], $s, $e, @F[4,5,6,7,9,10,11]);'  {output}.tmp2 | sort -k3,3 -k1gr > {output}
	rm {output}.tmp2
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
         bam="{genome}-MAP-{reads}-N.bam", bai="{genome}-MAP-{reads}-N.bam.bai"
    shell:
        """
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
    	samtools index {output.bam}
        """

# nanopore aligned by minimap for reads called nanopore.fastq.gz
rule Nanopore_minimap_bam2:
    input:
         fa="{genome}.fa", fastq="nanopore.fastq.gz"
    output:
         bam="{genome}-N.bam", bai="{genome}-N.bam.bai"
    shell:
        """
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
    	samtools index {output.bam}
        """

# nanopore aligned by minimap for reads called nanopore-sample.fastq.gz
rule Nanopore_minimap_bam_sample:
    input:
         fa="{genome}.fa", fastq="nanopore-sample.fastq.gz"
    output:
         bam="{genome}-NS.bam", bai="{genome}-NS.bam.bai"
    shell:
        """
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
    	samtools index {output.bam}
        """

# create coverage from bam file (for illumina / rnaseq coverage)
rule bam_to_bedgraph:
    input:
         bam="{aln}.bam"
    output:
         cov="{aln}.bedgraph"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bga -split > {output.cov}
        """

rule bedgraph_to_bw:
    input:
         bedgraph="{genome}-{aln}.bedgraph", sizes="{genome}.sizes"
    output:
         "{genome}-{aln}.bw"
    shell:
        """
	sort --buffer-size=1G -k1,1 -k2,2g {input.bedgraph} > {output}.tmp
        bedGraphToBigWig {output}.tmp {input.sizes} {output}
        rm {output}.tmp
        """


# nanopore aligned by minimap with paf format
rule Nanopore_minimap_paf:
    input:
         fa="{genome}.fa", fastq="{reads}-N.fastq.gz"
    output:
         paf="{genome}-MAP-{reads}-N.paf"
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# nanopore aligned by minimap with paf format for reads called nanopore.fastq.gz
rule Nanopore_minimap_paf2:
    input:
         fa="{genome}.fa", fastq="nanopore.fastq.gz"
    output:
         paf="{genome}-N.paf"
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# nanopore aligned by minimap with paf format for reads called nanopore-sample.fastq.gz
rule Nanopore_minimap_paf_sample:
    input:
         fa="{genome}.fa", fastq="nanopore-sample.fastq.gz"
    output:
         paf="{genome}-NS.paf"
    shell:
        """
        minimap2 -c -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} > {output.paf}
        """

# read connections between contigs
rule Nanopore_connections:
    input:
         fa="{genome}.fa", view="{genome}-{reads}.paf.view"
    output:
         "{genome}-{reads}.pairs"
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

# atom sequences in reads
# {name} is typically {genome}-{reads}
# do filtering on paf.view before running this
# column $F[2] is read name, $F[6] atom name
rule atom_sequences:
    input:
         "{name}.paf.view"
    output:
         "{name}.atoms"
    shell:
        """
        sort -k3,3 -k5,5g {input} > {output}.tmp1
        perl -lane 'die unless @F==11; if(!defined $o || $o ne $F[2]) {{ printf "\n%s", $F[2]; $o=$F[2]; }} printf " %s%s", $F[6], $F[1]; END {{ print ""; }} ' {output}.tmp1 | grep . > {output}
	rm {output}.tmp[1]
        """


# clipped nanopore reads
rule Nanopore_clipped_left:
    input:
         fa="{genome}.fa", view="{genome}-{reads}.paf.view"
    output:
         "{genome}-{reads}.clipL{size}.bed"
    shell:
        """
        perl -lane 'die unless @F==11; if($F[1] eq "+") {{ $gap=$F[4];}} else {{ $gap=$F[3]-$F[5]; }} if($gap >= {wildcards.size}) {{ print join("\t", $F[6], $F[8], $F[8]+1, $F[2], $gap, $F[1]); }} ' {input.view} | sort -k1,1 -k2,2g > {output}
        """

# clipped nanopore reads
rule Nanopore_clipped_right:
    input:
         fa="{genome}.fa", view="{genome}-{reads}.paf.view"
    output:
         "{genome}-{reads}.clipR{size}.bed"
    shell:
        """
        perl -lane 'die unless @F==11; if($F[1] eq "-") {{ $gap=$F[4];}} else {{ $gap=$F[3]-$F[5]; }} if($gap >= {wildcards.size}) {{ print join("\t", $F[6], $F[9], $F[9]+1, $F[2], $gap, $F[1]); }} ' {input.view} | sort -k1,1 -k2,2g > {output}
        """

# bed to bedgraph
rule bed_to_bedgraph:
    input:
         bed="{genome}-{aln}.bed", sizes="{genome}.sizes"
    output:
         cov="{genome}-{aln}.bedgraph"
    shell:
        """
        bedtools genomecov -i {input.bed} -bga -split -g {input.sizes} > {output.cov}
        """

# cluster by cd hit
# sequence names should be unique in first 19 chars
rule cdhit_cluster:
    input:
         "{name}.fa"
    output:
         fa="{name}-cl{id,[0-9\.]+}.fa",
         clstr="{name}-cl{id}.fa.clstr",
    shell:
        """
        cdhit-est -i {input} -o {output.fa} -c {wildcards.id}
        """

# compute cluster sizes from cdhit output
# the output will contain size of cluster for each cluster center
# sequence names should be unique in first 19 chars
rule cdhit_sizes:
    input:
         clstr="{name}.fa.clstr",
    output:
         sizes="{name}.fa.clsizes"
    shell:
        """
        perl -lane 'if(/^>/ && defined $n) {{ die unless defined $s; print "$s $n"; $n=0; $s=undef; }} elsif(/\*$/) {{ die unless />(.*)\s+\*$/; $s=$1; $s=~s/\.\.\.$//; $n++; }} else {{ $n++}} END {{  die unless defined $s; print "$s $n"; }}' < {input.clstr} > {output.sizes}
        """


rule cdhit_sizefilter:
    input:
         fa="{name}.fa",
	 sizes="{name}.fa.clsizes"
    output:
         fa="{name}-min{size,[0-9]+}.fa", sizes="{name}-min{size}.fa.clsizes",
    shell:
        """
        # get big cluster reps
        perl -lane 'die unless @F==2; if($F[1]>={wildcards.size}) {{ print }}' {input.sizes} > {output.sizes}
        # format for grep
	perl -lane 'print ">$F[0]"' < {output.sizes} > {output.fa}.tmp1
        # get names from fa for grep
        grep ">" {input.fa} > {output.fa}.tmp2
        # find selected names
        grep -F -f {output.fa}.tmp1 {output.fa}.tmp2 | perl -lne 's/>//; print' > {output.fa}.tmp3
        # get full fasta records
	faSomeRecords {input.fa} {output.fa}.tmp3 {output.fa}
        rm {output.fa}.tmp[123]
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

rule illumina_sample:
    input:
         fq1="{name}-I_1.fastq.gz", fq2="{name}-I_2.fastq.gz"
    output:
         fq1="{name}-sample-{frac}-I_1.fastq.gz", fq2="{name}-sample-{frac}-I_2.fastq.gz"
    params:
        list="{name}-sample-{frac}-I.list"
    shell:
        """
        zcat {input.fq1} | perl -lane 'BEGIN {{ srand({SAMPLE_SEED}); }} if($.%4==1 && rand(1)<{wildcards.frac}) {{ print $F[0]; }}' > {params.list}
        zcat {input.fq1} | grep -F -f {params.list} -A 3 - | grep -v '^--$' | gzip -c > {output.fq1}
        zcat {input.fq2} | grep -F -f {params.list} -A 3 - | grep -v '^--$' | gzip -c > {output.fq2}
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
        psl="{name}-self.psl", tab="{name}-self.psl.tab"
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

rule last_split_genome_aln:
    input:
        fa1="{name1}.fa", fa2="{name2}.fa"
    output:
        psl="{name1}-LASTSPLIT-{name2}.psl", tab="{name1}-LASTSPLIT-{name2}.psl.tab"
    shell:
        """
        lastdb {output.psl}-tmp {input.fa1}
        lastal {output.psl}-tmp {input.fa2} -E1e-10 | last-split > {output.psl}.maf
        maf-convert psl {output.psl}.maf > {output.psl}
        maf-convert tab {output.psl}.maf > {output.tab}
        rm {output.psl}.maf {output.psl}-tmp.*
        """


# align reads via last with last-split (starting from reads in fasta)
rule nano_last:
    input:
        ref="{genome}.fa", reads="{reads}.fa"
    output:
        psl="{genome}-MAP-{reads}.psl"
    shell:
        """
        lastdb -uNEAR {output.psl}-tmp {input.ref}
        lastal {output.psl}-tmp {input.reads} -E1e-10 | last-split | maf-convert psl  > {output.psl}
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

rule psl_view:
    input:
        "{name}.psl"
    output:
        view="{name}.psl.view", view2="{name}.psl.view2"
    shell:
        """
        perl -lane 'next unless /^[0-9]/; $g=$F[0]+$F[2]; $b=$F[1]+$F[3]+$F[5]+$F[7]; $s=sprintf "%.1f", 100*$g/($g+$b); print join("\t", @F[0..16], $s)' {input} > {output.view}
        perl -lane 'print join("\t", @F[0,8..17])' {output.view} > {output.view2}
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

rule freebayes:
    input:
        fa="{name}.fa",
        bam="{name}-I.bam"
    output: "{name}-I.vcf"
    shell: """
        freebayes                           \
            --bam {input.bam}               \
            --fasta-reference {input.fa}   \
            --vcf {output}
        # --min-mapping-quality 50 --min-base-quality 20 --min-alternate-count 3
        # --min-alternate-fraction 0.5 --min-coverage 3 --pooled-continuous
    """

rule freebayes2:
    input:
        fa="{name}.fa",
        bam="{name}-MAP-{reads}-I.bam"
    output: "{name}-MAP-{reads}-I.vcf"
    shell: """
        freebayes                           \
            --bam {input.bam}               \
            --fasta-reference {input.fa}   \
            --vcf {output}
        # --min-mapping-quality 50 --min-base-quality 20 --min-alternate-count 3
        # --min-alternate-fraction 0.5 --min-coverage 3 --pooled-continuous
    """


rule filter_vcf:
    input: "{name}-I.vcf",
    output: "{name}-If.vcf"
    shell: """
        bcftools view \
          --trim-alt-alleles \
          --exclude 'GT="ref" || QUAL<5' \
          {input} > {output}
    """

# polishing by freebayes
# take variants with genotype 1/1 or 1/2 and apply them
rule vcf_apply_alt_only:
    input:
      vcf="{name}-If.vcf",
      fa="{name}.fa",
    output:
      vcf="{name}-If1.vcf.gz",
      fa="{name}F.fa"
    shell: """
        bcftools view --include 'GT="1/1" || GT="1/2"' {input.vcf} | bgzip -c > {output.vcf}
        bcftools index -t {output.vcf}
	bcftools consensus -H 1 -f {input.fa} {output.vcf} > {output.fa}
    """

# create two haplotypes
#  preferably on sequence polished by freebayes (suffix F.fa)
rule hapcut2:
    input:
        bam="{name}-N.bam",
        vcf="{name}-If.vcf",
        fa="{name}.fa"
    log:
        hair="{name}H.hairlog",
        cut="{name}H.cutlog",
    output:
        bl="{name}H.blocks",
        vcfgz="{name}H.vcf.gz",
        tbi="{name}H.vcf.gz.tbi"
    params: "{name}H"
    shell: """
        extractHAIRS \
           --bam {input.bam} \
           --VCF {input.vcf} \
           --ont 1 \
           --indels 1 \
           --triallelic 1 \
           --ref {input.fa} \
           --out {params}.fragment 2> {log.hair}
        HAPCUT2 \
            --fragments {params}.fragment \
            --VCF {input.vcf} \
            --output {output.bl} \
            --outvcf 1 2> {log.cut}
        mv {output.bl}.phased.VCF {params}.vcf
        bgzip {params}.vcf
        bcftools index -t {output.vcfgz}  
    """

rule haplotyped_fasta:
    input:
        fa="{name}.fa",
        vcf="{name}H.vcf.gz",
        tbi="{name}H.vcf.gz.tbi"
    output: "{name}H{haplotype,[12]}.fa"
    shell: """
        bcftools consensus -H {wildcards.haplotype} -f {input.fa} {input.vcf} > {output}
    """

rule concat_haplotyped_fasta:
   input:
     fa1="{name}H1.fa",
     fa2="{name}H2.fa"
   output:
     "{name}H.fa"
   shell:
     """
     perl -lne 's/>(.*)$/>$1_h1/; print;' {input.fa1} > {output}.tmp1 
     perl -lne 's/>(.*)$/>$1_h2/; print;' {input.fa2} > {output}.tmp2
     cat {output}.tmp1 {output}.tmp2 > {output}
     rm {output}.tmp1 {output}.tmp2
     """

rule kmer_bedgraph:
   input: "{name}.fa"
   output: "{name}-kmers{K,[0-9]+}.bedgraph"
   shell:
      """
      {SCRIPT_PATH}/find_kmers.pl 21 < {input} | sort --buffer-size=1G | uniq -c > {output}.tmp1
      {SCRIPT_PATH}/find_kmers.pl -b 21 < {input} | sort --buffer-size=1G -k4 > {output}.tmp2
      join -1 2 -2 4 {output}.tmp1 {output}.tmp2 | perl -lane 'print join("\t", @F[2,3,4,1])' | sort --buffer-size=1G -k1,1 -k2,2g > {output}
      rm {output}.tmp1 {output}.tmp2
      """

# create sliding windows, keep only full-length ones
# (except for full sequence of too short)
rule windows:
    input: "{name}.sizes"
    output: "{name}-windows{size}bp-{slide}bp.bed"
    shell:
      """      
      bedtools makewindows -w {wildcards.size} -s {wildcards.slide} -g {input} | perl -lane 'print if $F[1]==0 || $F[2]-$F[1]=={wildcards.size}' > {output}
      """

# take median in each windw, slow
rule smooth_bedgraph:
     input:
       bg="{genome}-{name}.bedgraph",
       bed="{genome}-windows{size}bp-1bp.bed"
     output: "{genome}-{name}-smooth{size}bp.bedgraph"
     shell:
       """
       bedtools map -b {input.bg} -a {input.bed} -c 4 -o median > {output}.tmp
       perl -lane '$h=int(($F[2]+$F[1])/2); print join("\t", $F[0], $h, $h+1, $F[3]);' {output}.tmp > {output}
       rm {output}.tmp
       """

# compute overall coverage median
rule bedgraph_median:
     input: "{name}.bedgraph"
     output: "{name}.bedgraph-median"
     shell:
       """
       R --vanilla --slave -e 'a=read.table("{input}"); cat(median(rep(a$V4, times=a$V3-a$V2)));' > {output}
       """

# compute coverage normalized by median
rule bedgraph_norm:
     input:
       bg="{name}.bedgraph",
       med="{name}.bedgraph-median"
     output: "{name}-norm.bedgraph"
     shell:
       """
       perl -lane 'BEGIN {{ $X = `cat {input.med}`; }} $F[3] = sprintf("%.5f", $F[3]/$X); print join("\\t", @F); ' {input.bg} > {output}
       """
