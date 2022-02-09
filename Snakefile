"""``snakemake`` file that runs analysis."""


import itertools
import os


configfile: 'config.yaml'


wildcard_constraints:
    accession='SRR\d+|all\-accessions'


rule all:
    input:
        'results/alignment_counts.csv',
        'results/variant_analysis.csv',
        'results/variant_analysis_to_outgroup.csv',

rule get_sra:
    """Download ``*.sra`` files."""
    output: sra=protected("results/sra_downloads/{accession}.sra")
    conda: 'environment.yml'
    shell: "prefetch {wildcards.accession} -o {output.sra}"

rule sra_fastq:
    """Extract FASTQ from SRA file."""
    input: sra=rules.get_sra.output.sra
    output: fastq="results/fastq_from_sra/{accession}.fastq.gz"
    conda: 'environment.yml'
    shell:
        """
        fastq-dump \
            --skip-technical \
            --readids \
            --read-filter pass \
            --dumpbase \
            --split-spot \
            --clip \
            --stdout \
            {input.sra} | \
        gzip > {output.fastq}
        """

rule preprocess_fastq:
    """Pre-process FASTQ with ``fastp``."""
    input: fastq=rules.sra_fastq.output.fastq
    output:
        fastq="results/preprocessed_fastqs/{accession}.fastq.gz",
        json="results/preprocessed_fastqs/{accession}.json",
        html="results/preprocessed_fastqs/{accession}.html",
    conda: 'environment.yml'
    shell:
        """
        fastp \
            --thread 1 \
            -i {input.fastq} \
            -q 15 \
            -u 40 \
            -l 25 \
            --trim_poly_g \
            --trim_poly_x \
            -o {output.fastq} \
            -j {output.json} \
            -h {output.html}
        """

rule get_refgenome:
    """Get reference genome."""
    output: fasta="results/refgenome/refgenome_untrimmed.fa"
    params: ftp=config['refgenome']
    conda: 'environment.yml'
    shell: "wget -O - {params.ftp} | gunzip -c > {output.fasta}"

rule refgenome_trim3_polyA:
    """Trim 3' polyA from reference genome."""
    input: fasta=rules.get_refgenome.output.fasta
    output: fasta="results/refgenome/refgenome_trimmed.fa"
    conda: 'environment.yml'
    script: 'scripts/trim3_polyA.py'

rule minimap2_genome:
    """Build ``minimap2`` reference genome."""
    input: fasta=rules.refgenome_trim3_polyA.output.fasta
    output: mmi="results/refgenome/refgenome_trimmed.mmi"
    conda: 'environment.yml'
    shell: "minimap2 -x {config[minimap2_setting]} -d {output.mmi} {input.fasta}"

rule align_fastq:
    """Align FASTQ to reference genome."""
    input:
        fastq=rules.preprocess_fastq.output.fastq,
        mmi=rules.minimap2_genome.output.mmi,
    output:
        sam=temp("results/alignments/{accession}.sam"),
        unsorted_bam=temp("results/alignments/{accession}.bam"),
        bam="results/alignments/{accession}_sorted.bam",
    conda: 'environment.yml'
    shell:
        """
        minimap2 \
            -a \
            -MD \
            -c \
            -eqx \
            -x {config[minimap2_setting]} \
            --sam-hit-only \
            --secondary=no \
            {input.mmi} \
            {input.fastq} \
            > {output.sam}
        samtools view -b -F 4 -o {output.unsorted_bam} {output.sam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        """

rule merge_alignments:
    """Aggregate alignments BAMs"""
    input: bams=expand(rules.align_fastq.output.bam, accession=config['accessions'])
    output: bam='results/alignments/all-accessions_sorted.bam'
    conda: 'environment.yml'
    shell: "samtools merge {output.bam} {input.bams}"

rule high_ident_bam:
    """Get BAM with just high-identity reads."""
    input: bam="results/alignments/{accession}_sorted.bam"
    output: bam="results/alignments/{accession}_high_ident_sorted.bam"
    params: **config['high_ident_thresholds']
    conda: 'environment.yml'
    script: 'scripts/high_ident_bam.py'

rule alignment_counts:
    """Aggregate alignment counts across accessions."""
    input:
        bams=expand(rules.align_fastq.output.bam, accession=config['accessions']),
        high_ident_bams=expand(rules.high_ident_bam.output.bam,
                               accession=config['accessions']),
    output: csv='results/alignment_counts.csv'
    conda: 'environment.yml'
    script: 'scripts/alignment_counts.py'

rule consensus_seq:
    """Create consensus sequence using ``ivar``."""
    input: bam="results/alignments/{bambase}_sorted.bam"
    output: fasta="results/consensus_seqs/{bambase}.fa"
    params: prefix=lambda wc, output: os.path.splitext(output.fasta)[0],
    conda: 'environment.yml'
    shell:
        """
        samtools mpileup -aa -A -d 0 -Q 20 {input.bam} | \
        ivar consensus \
            -p {params.prefix} \
            -t 0.5 \
            -m 1
        """

rule ivar_variants:
    """Call variants using ``ivar``."""
    input:
        bam="results/alignments/{bambase}_sorted.bam",
        ref=rules.refgenome_trim3_polyA.output.fasta,
    output: fasta="results/variants/{bambase}.tsv"
    params: prefix=lambda wc, output: os.path.splitext(output.fasta)[0],
    conda: 'environment.yml'
    shell:
        """
        samtools mpileup -aa -A -d 0 -B -Q 20 {input.bam} | \
        ivar variants \
            -p {params.prefix} \
            -t 0.5 \
            -r {input.ref}
        """

rule outgroup_fasta:
    """Get FASTA for outgroup."""
    output: fasta="results/outgroup/{outgroup}.fa"
    params: acc=lambda wc: config['outgroups'][wc.outgroup]
    conda: 'environment.yml'
    shell:
        """
        efetch \
            -format fasta \
            -db nuccore \
            -id {params.acc} \
            > {output.fasta}
        """

rule outgroup_alignment:
    """Align outgroup to reference."""
    input:
        ref=rules.refgenome_trim3_polyA.output.fasta,
        outgroups=expand(rules.outgroup_fasta.output.fasta,
                         outgroup=config['outgroups'])
    output:
        concat_fasta=temp('results/outgroup/to_align.fa'),
        alignment='results/outgroup/alignment.fa',
    conda: 'environment.yml'
    shell:
        # insert newline between FASTA files when concatenating:
        # https://stackoverflow.com/a/25030513/4191652
        """
        awk 1 {input} > {output.concat_fasta}
        mafft {output.concat_fasta} > {output.alignment}
        """

rule outgroup_map:
    """Map sites in viral genome to outgroup identities."""
    input: alignment=rules.outgroup_alignment.output.alignment
    output: site_map='results/outgroup/site_map.csv'
    conda: 'environment.yml'
    script: 'scripts/outgroup_map.py'

rule analyze_variants:
    """Analyze variants present in sequences."""
    input:
        site_map=rules.outgroup_map.output.site_map,
        variants=expand("results/variants/{accession}{high_ident}.tsv",
                        accession=config['accessions'] + ['all-accessions'],
                        high_ident=['', '_high_ident']),
    output:
        all_csv='results/variant_analysis.csv',
        to_outgroup_csv='results/variant_analysis_to_outgroup.csv',
    params:
        descriptors=[{'accession': accession, 'identity': ident}
                     for (accession, ident) in
                     itertools.product(config['accessions'] + ['all-accession'],
                                       ['any', 'high'])
                     ]
    conda: 'environment.yml'
    script: 'scripts/analyze_variants.py'
