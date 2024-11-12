MINIMAP_OPT = "-t 4"

import os
if os.path.exists("nanopore.yaml"):
   config_path = "nanopore.yaml"
else:
  config_path = f"{workflow.basedir}/nanopore.yaml"
  assert os.path.exists(config_path)

configfile: config_path


# basecall a single pod5 file using fast basecaller
rule pod5ToBamFast:
     input: "{name}.pod5"
     output: "{name}-fast-unaligned.bam"
     shell:
       """
       {config[dorado_path]} basecaller --device {config[dorado_device]} '{config[dorado_model_fast_path]}' {input} > {output}
       """ 

# basecall whole reads folder using slow basecaller
rule pod5ToBamSlow:
     input: directory("reads")
     output: "reads-unaligned.bam"
     shell:
       """
       {config[dorado_path]} basecaller --device {config[dorado_device]} '{config[dorado_model_slow_path]}' {input}/*/*/pod5 > {output}
       """ 



# convert bam to fastq
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

# nanopore aligned to ref.fa by minimap with bam format
rule Nanopore_minimap_bam_ref:
    input:
         fa="ref.fa", fastq="{reads}.fastq.gz"
    output:
         bam="{reads}-mapped.bam", bai="{reads}-mapped.bam.bai"
    shell:
        """
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.fa} {input.fastq} | samtools view -S -b - | samtools sort - -o {output.bam}
	samtools index {output.bam}
        """


# nanoplot on nanopore reads in fastq
rule nanoplotFastq:
     input: "{name}.fastq.gz"
     output:
        dir=directory("{name}.fastq.len{size}.plot"),
        txt="{name}.fastq.len{size}.plot.txt"
     shell:
        """
        NanoPlot --fastq {input} --maxlength {wildcards.size} -o {output.dir}
        cp {output.dir}/NanoStats.txt {output.txt}
        """

# nanoplot on nanopore reads in bam
rule nanoplotBam:
     input: "{name}-mapped.bam"
     output:
        dir=directory("{name}-mapped.bam.len{size}.plot"),
        txt="{name}-mapped.bam.len{size}.plot.txt"
     shell:
        """
        samtools view -b -F 2048 {input} > {output.txt}.tmp.bam
        NanoPlot --tsv_stats --bam {output.txt}.tmp.bam --maxlength {wildcards.size} -o {output.dir}
        cp {output.dir}/NanoStats.txt {output.txt}
        rm {output.txt}.tmp.bam
        rm -f {output.txt}.tmp.bam.bai
        """

# nanoplot on nanopore reads in unaligned bam
rule nanoplotUBam:
     input: "{name}-unaligned.bam"
     output:
        dir=directory("{name}-unaligned.bam.len{size}.plot"),
        txt="{name}-unaligned.bam.len{size}.plot.txt"
     shell:
        """
        NanoPlot --tsv_stats --ubam {input} --maxlength {wildcards.size} -o {output.dir}
        cp {output.dir}/NanoStats.txt {output.txt}
        """

# overall stats with length specified in config
# length influences only the plots, not stats
rule new_stats:
  input: "{name}.len" + str(config['nanoplot_maxlen']) + ".new_stats"
  output: "{name}.new_stats"
  shell:
    """
    cp {input} {output}
    """

nanoplot_columns = [
  ('reads','number_of_reads'),
  ('total_bp','number_of_bases'),
  ('median_len', 'median_read_length'),
  ('mean_len', 'mean_read_length'),
  ('n50_len', 'n50'),
  ('med_id', 'median_identity'),
  ('med_qual', 'median_qual')
]

rule new_stats_len:
  input:
    stu="{name}-unaligned.bam.len{size}.plot.txt",
    stm="{name}-mapped.bam.len{size}.plot.txt"
  output:
    "{name}.len{size}.new_stats"
  run:
    import csv
    from collections import defaultdict
    with open(output[0], "w") as out:
      d = dict()
      out.write(" ".join(['      '] + [x[0] for x in nanoplot_columns]) + "\n")
      for (which, filename) in [('all   ', input.stu), ('mapped', input.stm)]:
        d[which] = defaultdict(lambda: "----")
        with open(filename, "r") as infile:
           reader = csv.reader(infile, delimiter='\t', quotechar='"')
           for row in reader:
              d[which][row[0]] = row[1]
        items = [which] + [d[which][x[1]] for x in nanoplot_columns]
        out.write("  ".join(items) + "\n")
    

  


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
