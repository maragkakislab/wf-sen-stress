# Define parameters. Will use defaults if cannot get them from config
NANOPOLISH_DIR = config.get('NANOPOLISH_DIR', 'nanopolish')
NANOPOLISH_RES = config.get('NANOPOLISH_RES', NANOPOLISH_DIR + '/results')


localrules: nanopolish_all, nanopolish_all_ebola


rule nanopolish_all:
    input:
        expand(
                NANOPOLISH_RES + "/{s}/polya.filtered.tab",
                s = [s.name for s in samples.values() if s.is_direct_rna()]),


rule nanopolish_all_ebola:
    input:
        expand(
                NANOPOLISH_RES + "/{s}/ebola.polya.filtered.tab",
                s = [s.name for s in samples.values() if s.is_direct_rna()]),
        expand(
                NANOPOLISH_RES + "/{s}/ebola.polya.filtered.fft.pdf",
                s = [s.name for s in samples.values() if s.is_direct_rna()]),


# nanopolish_index creates an index mapping from basecalled reads in FASTQ to
# the original fast5s containing each read.
rule nanopolish_index:
    input:
        origin = lambda ws: EXP_DIR + '/' + samples[ws.s].parent_exp
                             + '/origin.txt',
        seq_sum = lambda ws: EXP_DIR + '/' + samples[ws.s].parent_exp
                             + '/guppy/sequencing_summary.txt',
        fastq = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz",
    output:
        SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index",
        SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.fai",
        SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.gzi",
        SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.readdb",
    params:
        fast5_dir = lambda wilds, input: os.path.dirname(input.origin)
                                         + "/runs",
    threads: 1
    resources:
        mem_mb=5*1024,
        runtime=1*24*60
    envmodules:
        "nanopolish/0.14.0.conda"
    shell:
        """
        nanopolish index \
                --directory {params.fast5_dir} \
                --sequencing-summary {input.seq_sum} \
                {input.fastq}
        """


# nanopolish_polya runs nanopolish to calculate the length of the poly(A) tail
# for each read
rule nanopolish_polya:
    input:
        fastq = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz",
        bam = SAMPLES_DIR + "/{s}/align/reads.toTranscriptome.sorted.bam",
        genome = DATA_DIR + '/' + ASSEMBLY + "/transcripts.fa",
        index = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index",
        index_fai = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.fai",
        index_gzi = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.gzi",
        index_readdb = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.readdb",
    output:
        NANOPOLISH_RES + "/{s}/polya.tab"
    threads: 20
    resources:
        runtime=4*24*60
    conda:
        "../envs/bam.yml"
    envmodules:
        "nanopolish/0.14.0.conda",
    shell:
        """
        nanopolish polya \
            --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --threads {threads} \
            | table-paste-col.py \
                --table - \
                --col-name sample \
                --col-val {wildcards.s} \
            > {output}
        """


# nanopolish_polya runs nanopolish to calculate the length of the poly(A) tail
# for each read
rule nanopolish_polya_ebola:
    input:
        fastq = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz",
        bam = SAMPLES_DIR + "/{s}/align/reads.toEbolaGenome.sorted.bam",
        genome = DATA_DIR + "/ebola/genome.fa",
        index = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index",
        index_fai = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.fai",
        index_gzi = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.gzi",
        index_readdb = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.readdb",
    output:
        NANOPOLISH_RES + "/{s}/ebola.polya.tab"
    threads: 20
    resources:
        runtime=4*24*60
    conda:
        "../envs/bam.yml"
    envmodules:
        "nanopolish/0.14.0.conda",
    shell:
        """
        nanopolish polya \
            --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --threads {threads} \
            | table-paste-col.py \
                --table - \
                --col-name sample \
                --col-val {wildcards.s} \
            > {output}
        """


# nanopolish_eventalign runs the eventalign command of nanopolish.
rule nanopolish_eventalign:
    input:
        fastq = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz",
        bam = SAMPLES_DIR + "/{s}/align/reads.toTranscriptome.sorted.bam",
        genome = DATA_DIR + '/' + ASSEMBLY + "/transcripts.fa",
        index = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index",
        index_fai = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.fai",
        index_gzi = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.gzi",
        index_readdb = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.readdb",
    output:
       summary = NANOPOLISH_RES + "/{s}/toTranscriptome.summary.txt",
       ea = temp(NANOPOLISH_RES + "/{s}/toTranscriptome.eventalign.txt")
    threads: 20
    resources:
        runtime=4*24*60
    envmodules:
        "nanopolish/0.14.0.conda",
    conda:
        "../envs/nanopolish.yml"
    shell:
        """
        nanopolish eventalign \
            --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --scale-events \
            --signal-index \
            --summary {output.summary}  \
            --threads {threads} \
            > {output.ea}
        """


# nanopolish_eventalign_ebola runs the eventalign command of nanopolish.
rule nanopolish_eventalign_ebola:
    input:
        fastq = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz",
        bam = SAMPLES_DIR + "/{s}/align/reads.toEbolaGenome.sorted.bam",
        genome = DATA_DIR + "/ebola/genome.fa",
        index = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index",
        index_fai = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.fai",
        index_gzi = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.gzi",
        index_readdb = SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz.index.readdb",
    output:
       summary = NANOPOLISH_RES + "/{s}/toEbolaGenome.summary.txt",
       ea = temp(NANOPOLISH_RES + "/{s}/toEbolaGenome.eventalign.txt")
    threads: 20
    resources:
        runtime=4*24*60
    envmodules:
        "nanopolish/0.14.0.conda",
    conda:
        "../envs/nanopolish.yml"
    shell:
        """
        nanopolish eventalign \
            --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --scale-events \
            --signal-index \
            --summary {output.summary}  \
            --threads {threads} \
            > {output.ea}
        """


# nanopolish_filter_passes filters the nanopolish poly(A) output to only keep
# reads with enough evidence.
rule nanopolish_filter_passes:
    input:
        "{filepath}polya.tab"
    output:
        "{filepath}polya.filtered.tab"
    shell:
        """
        cat {input} \
            | grep -P '(qc_tag|PASS)' \
            > {output}
        """


# run_fft_for_polya runs Fourier Transform for the distribution of poly(A)
# tail lengths.
rule run_fft_for_polya:
    input:
        polya = ANALYSIS_DIR + "/polya/results/{s}/ebola.polya.filtered.tab",
        bed = DATA_DIR + "/ebola/genes.mRNA.bed",
        bam = SAMPLES_DIR + "/{s}/align/reads.toEbolaGenome.sorted.bam",
    output:
        NANOPOLISH_RES + "/{s}/ebola.polya.filtered.fft.pdf"
    shell:
        """
        python3 workflow/scripts/nanopolish/plot_polya_dist_for_feats.py \
                --bed {input.bed} \
                --bam {input.bam} \
                {input.polya} \
                {output}
        """


# nanopolish_repack_vbz_to_gzip changes the compression of the fast5 files
# from vbz to gzip. vbz is a new compression algorithm used by ONT but is
# creating problems with nanopolish polya.
rule nanopolish_repack_vbz_to_gzip:
    input:
        origin = EXP_DIR + '/{e}/origin.txt',
        seq_sum = EXP_DIR + '/{e}/guppy/sequencing_summary.txt',
    output:
        odir = temp(directory(EXP_DIR + '/{e}/repack')),
        sentinel = EXP_DIR + '/{e}/repack/repacked.txt',
        seq_sum = EXP_DIR + '/{e}/repack/sequencing_summary.txt',
    params:
        fast5_dir = lambda wilds, input: os.path.dirname(input.origin),
    threads: 40
    resources:
        mem_mb=10*1024,
        runtime=1*24*60
    envmodules:
        'ont-fast5-api/4.0.0',
    shell:
        """
        compress_fast5 \
                --input_path {params.fast5_dir}/runs/ \
                --save_path {params.fast5_dir}/repack/ \
                --compression gzip \
                --threads {threads} \
                --recursive

        cp {input.seq_sum} {output.seq_sum}

        touch {output.sentinel}
        """
