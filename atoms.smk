# atom sequences for genomes:
# {genome}.fa, {atoms}.fa
# {genome}-MAPF-{atoms}.atomsBoth
# {genome}-MAPF-{atoms}.atomGapsBoth

# atoms sequences for reads in nanopore.fastq.gz
# wrt {atoms}.fa
# {atoms}-Naf.atomGapsBoth
# {atoms}-Naf.atomsBoth

# folder with assembled sequences for reads containing individual atoms
# or atom combinations
# create file {name}-atomassembly.txt
#   each line has some id, tab, regular expression to be applied to rows of atomsBoth
# link {atoms}-Naf.atomsBoth as {name}-atomassembly.atoms
# link {atoms}.fa as {name}-atomassembly.atoms.fa
# first create {name}-atomassembly-dir
# then create {name}-atomassembly-dir/done


# align nanopore reads to atoms
#   special settings to increase telo alignments
#   -f filter out top FLOAT fraction of repetitive minimizers 0.0002->0
#   -w minimizer window size 10->5
rule atom_aln:
     input: fa="{name}.fa", fq="nanopore.fastq.gz"
     output: paf="{name}-Na.paf"
     shell:
       """
       minimap2 -c -x map-ont -w 5 -f 0 --secondary=no -t 4 {input.fa} {input.fq} > {output.paf}
       """

# align assembly to atoms
#   similarly as reads
rule assembly_aln:
     input: atoms="{atoms}.fa", fa="{name}.fa"
     output: paf="{name}-MAP-{atoms}.paf"
     shell:
       """
       minimap2 -c -x map-ont -w 5 -f 0 --secondary=no -t 4 {input.atoms} {input.fa} > {output.paf}
       """

# create a more readable version of a paf file
rule minimap_view:
    input:
        "{aln}.paf"
    output:
        "{aln}.paf.view"
    shell:
        """
        perl -lane '$id=sprintf("%.1f", $F[9]*100/$F[10]); print join("\\t", @F[9,4,0..3,5..8],$id)' {input} > {output}
        """
# convert paf.view to bed
# from paf.view alignment file in which special sequences of interest
# were used as db and genome as query
# create a bed file for the genome containing name of special sequence
# together with % coverage and  %id; number of matches are used as score
rule paf_view_bed:
    input: "{name}.paf.view"
    output: "{name}.paf.bed"
    shell:
      """
      perl -lane '$cov=sprintf "%.1f", ($F[9]-$F[8])*100/$F[7]; $n=join("_", $F[6], "cov".$cov, "id".$F[10]); print join("\t", @F[2,4,5],$n,@F[0,1])' {input} | sort -k1,1 -k2g > {output}
      """

rule atom_filter:
     input: "{name}-Na.paf.view"
     output: "{name}-Naf.paf.view"
     shell:
       """
       perl -lane 'print if $F[0]>=0.5*$F[7]' {input} > {output}
       """

rule assembly_filter:
     input: "{name}-MAP-{atoms}.paf.view"
     output: "{name}-MAPF-{atoms}.paf.view"
     shell:
       """
       perl -lane 'print if $F[0]>=0.5*$F[7]' {input} > {output}
       """

#######

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
        perl -lane 'die unless @F==11; if(!defined $o || $o ne $F[2]) {{ printf "\\n%s", $F[2]; $o=$F[2]; }} printf " %s%s", $F[6], $F[1]; END {{ print ""; }} ' {output}.tmp1 | grep . > {output}
        rm {output}.tmp[1]
        """

# gaps between atoms
rule atom_gaps:
    input:
         "{name}.paf.view"
    output:
         "{name}.atomGaps"
    shell:
        """
        sort -k3,3 -k5,5g {input} > {output}.tmp1
        perl -lane 'die unless @F==11; if(!defined $o || $o ne $F[2]) {{ if(defined $g) {{ printf " %d\\n", $g; }} printf "%s", $F[2]; $o=$F[2]; $e=0; }} printf " %d %s%s", $F[4]-$e, $F[6], $F[1]; $e=$F[5]; $g=$F[3]-$F[5]; END {{ print " $g"; }} ' {output}.tmp1  > {output}
        rm {output}.tmp[1]
        """

# run-length encoding of atom sequences
# for easier analysis of tandem repeats
rule reduce_atom_sequences:
    input:
         "{name}.atoms"
    output:
         "{name}.atomsR"
    shell:
        """
        perl -lane 'my @G; $o=""; push @F, "X"; foreach $f (@F) {{ if($f eq $o) {{ $n++; }} else {{ if($n>1) {{ $o.="X".$n; }} push @G, $o; $n=1; $o=$f; }} }}; print join(" ", @G);' {input} > {output}
        """

# reverting atomGaps sequences can be used
# for normalizing them to one direction
rule revert_atomGaps_sequences:
    input:
         "{name}.atomGaps"
    output:
         "{name}.atomGapsRC"
    shell:
        """
        perl -lane '$n=shift @F; $n.="RC"; @F = reverse(@F); foreach $f (@F) {{if($f=~/\\+$/) {{$f=~s/\\+$/-/ or die; }} elsif ($f=~/\\-$/) {{ $f=~s/\\-$/+/ or die; }} else {{ die unless $f=~/^[0-9-]+$/; }} }} print join(" ", $n, @F); ' {input} > {output}
        """


# reverting atom sequences can be used
# for normalizing them to one direction
rule revert_atom_sequences:
    input:
         "{name}.atoms"
    output:
         "{name}.atomsRC"
    shell:
        """
        perl -lane '$n=shift @F; $n.="RC"; @F = reverse(@F); foreach $f (@F) {{if($f=~/\\+$/) {{$f=~s/\\+$/-/ or die; }} else {{ $f=~s/\\-$/+/ or die; }} }} print join(" ", $n, @F); ' {input} > {output}
        """

# combining both directions of atoms
rule combine_atom_sequences:
    input:
         "{name}.atoms", "{name}.atomsRC"
    output:
         "{name}.atomsBoth"
    shell:
        """
	cat {input} | sort > {output}
        """


# combining both directions of atom gaps
rule combine_atom_gaps:
    input:
         "{name}.atomGaps", "{name}.atomGapsRC"
    output:
         "{name}.atomGapsBoth"
    shell:
        """
	cat {input} | sort > {output}
        """


############

#use with suffix s (-atoms) or assembly (-atomassembly)
rule dir_prepare:
     input: "{name}-atom{suffix}.txt"
     output: directory("{name}-atom{suffix}-dir")
     run:
        import os
        dir=output[0]
        os.makedirs(dir, exist_ok=True)

        with open(input[0], 'r') as f:
            for line in f:
                parts = line.strip().split("\t")
                assert len(parts) >= 2, "Each line must contain at least 2 parts separated by tabs"
                (name, pattern) = parts[0], parts[1]
        
                output_dir = os.path.join(dir, f"{name}-subdir")
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, f"pattern.txt")
                with open(output_file, 'w') as out:
                     out.write(f"{pattern}\n")   


rule get_atom_seq:
     input: txt="{group}-atoms-dir/{name}/pattern.txt",
            atoms="{group}-atoms.atoms"
     output: atoms="{group}-atoms-dir/{name}/atoms.txt"
     shell:
       """
       grep -E -f {input.txt} -o {input.atoms} > {output.atoms} || true
       """

rule get_atom_prefixes:
     input: "{path}/atoms.txt"
     output: "{path}/atoms-prefixes.txt"
     shell:
       """
       perl -lane 'foreach $n (0..$#F) {{ print join(" ", $n+1, @F[0..$n]); }}' {input} > {output}
       """

rule get_atom_prefix_counts:
     input: "{path}/atoms-prefixes.txt"
     output: "{path}/atoms-prefixes-counts.txt"
     shell:
       """
       sort {input} | uniq -c | sort -k1,1gr -k2,2gr > {output}
       """


rule atomas_one_done:
     input: "{group}-atoms-dir/{name}/atoms-prefixes-counts.txt"
     output: "{group}-atoms-dir/{name}/done",
     shell:
        """
	touch {output}
	"""



#############
#
# 

def assembly_done_inputs(wildcards):
    dir = wildcards.group + "-atom" + wildcards.suffix + "-dir"
    IDS, = glob_wildcards(dir + "/{id}-subdir/")
    return expand(dir + "/{id}-subdir/done", id=IDS)
    

#use with suffix s (-atoms) or assembly (-atomassembly)
rule assembly_done_all:
     input: assembly_done_inputs
     output: "{group}-atom{suffix}-dir/done"
     shell:
        """
	echo {input}
	touch {output}
	"""

rule atom_get_read_list:
     input: txt="{group}-atomassembly-dir/{name}/pattern.txt",
            atoms="{group}-atomassembly.atoms"
     output: list="{group}-atomassembly-dir/{name}/nanopore.list"
     shell:
       """
       grep -E -f {input.txt} {input.atoms} | perl -lane '$F[0]=~s/RC$//; print $F[0]' > {output.list}
       """

rule atom_get_reads:
     input: fq="nanopore.fastq.gz",
            list="{group}-atomassembly-dir/{name}/nanopore.list"
     output: fq="{group}-atomassembly-dir/{name}/nanopore.fastq.gz",
     shell:
       """
       zcat {input.fq} | grep -F -f {input.list} -A 3 - | grep -v '^--$' | gzip -c > {output.fq}
       """

rule atom_flye:
     input: fq="{group}-atomassembly-dir/{name}/nanopore.fastq.gz"
     output: fa="{group}-atomassembly-dir/{name}/flye.fa"
     params: dir="{group}-atomassembly-dir/{name}"
     shell:
       """
       flye -t 4 --nano-raw {input.fq} --out-dir {params.dir}/flye
       mv {params.dir}/flye/assembly.fasta {output.fa}
       """

rule align_atoms_to_fa:
     input: fa="{group}-atomassembly-dir/{name}/flye.fa",
            atoms="{group}-atomassembly.atoms.fa"
     output: paf="{group}-atomassembly-dir/{name}/flye.paf"
     shell:
        """
        minimap2 -c -x map-ont -w 5 -f 0 --secondary=no -t 4 {input.atoms} {input.fa} > {output.paf}
        """

rule align_atoms_to_fa_bed:
     input: view="{group}-atomassembly-dir/{name}/flye.paf.view",
     output: bed="{group}-atomassembly-dir/{name}/flye.paf.bed",
     shell:
       """
       /opt/assembly-scripts/annot {output.bed}
       """

rule align_reads_to_fa_paf:
     input: fa="{group}-atomassembly-dir/{name}/flye.fa",
            fq="{group}-atomassembly-dir/{name}/nanopore.fastq.gz",
     output: paf="{group}-atomassembly-dir/{name}/flye-N.paf",
     	     view="{group}-atomassembly-dir/{name}/flye-N.paf.view"
     params: dir="{group}-atomassembly-dir/{name}"
     shell:
        """
	cd {params.dir}; /opt/assembly-scripts/assembly flye-N.paf.view
	"""

rule align_reads_to_fa_bam:
     input: fa="{group}-atomassembly-dir/{name}/flye.fa",
            fq="{group}-atomassembly-dir/{name}/nanopore.fastq.gz",
     output: bam="{group}-atomassembly-dir/{name}/flye-N.bam",
     	     bai="{group}-atomassembly-dir/{name}/flye-N.bam.bai"
     params: dir="{group}-atomassembly-dir/{name}"
     shell:
        """
	cd {params.dir}; /opt/assembly-scripts/assembly flye-N.bam.bai
	"""

rule atomassembly_one_done:
     input: bam="{group}-atomassembly-dir/{name}/flye-N.bam",
            bed="{group}-atomassembly-dir/{name}/flye.paf.bed",
            view="{group}-atomassembly-dir/{name}/flye.paf.view"
     output: "{group}-atomassembly-dir/{name}/done",
     shell:
        """
	touch {output}
	"""
