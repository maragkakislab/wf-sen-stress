DORADO_MODELS = config.get('DORADO_MODELS', {})
DORADO_DEFAULT_MODEL = config.get('DORADO_DEFAULT_MODEL', 'hac')

# # Rule dorado_fast5_to_pod5 converts fast5 files to pod5.
# rule dorado_fast5_to_pod5:
#     input:
#         RAW_DIR + "/{e}/origin.txt"
#     output:
#         directory(EXP_DIR + "/{e}/runs_pod5")
#     threads:
#         5 # NEVER set this higher on biowulf; crashes due to IO
#     resources:
#         mem_mb = 32*1024,
#         runtime = 8*24*60
#     conda:
#         "../envs/pod5.yml"
#     shell:
#         """
#         mkdir -p {output} 
#         pod5 convert fast5 \
#                 {RAW_DIR}/{wildcards.e}/runs/ \
#                 --recursive \
#                 --threads {threads} \
#                 --output {output}/
#         """


# dorado_barcode_options returns the additional options that are used in the
# dorado_basecaller to enable detection of barcodes.
def dorado_barcode_options(experiment):
    if experiment.is_barcoded:
        return '--kit-name ' + experiment.kit

    return ''


# dorado_unstranded_options returns the additional options that are used in
# the dorado_basecaller to disable detection of adaptors. Adaptors are used in
# subsequent steps (pychopper) to orient the reads.
def dorado_unstranded_options(experiment):
    if s.is_unstranded():
        return '--no-trim'

    return ''


# Rule dorado_basecall runs the dorado basecaller.
rule dorado_basecall:
    input:
        origin = RAW_DIR + "/pod5/{e}/origin.txt",
        file_dir = RAW_DIR + "/pod5/{e}/pod5"

    output:
        EXP_DIR + "/{e}/dorado/calls.bam"
    params:
        model = lambda wilds: DORADO_MODELS.get(wilds.e, DORADO_DEFAULT_MODEL),
        additional_opts = lambda wilds:
            " ".join([dorado_barcode_options(experiments[wilds.e]),
                      dorado_unstranded_options(experiments[wilds.e])]),
    threads:
        8
    resources:
        gpu = 2,
        gpu_model = "[gpuv100x|gpua100]",
        mem_mb = 64*1024,
        runtime = 1*24*60
    envmodules:
        "dorado/0.7.1"
    shell:
        """
        dorado basecaller \
            --recursive \
            {params.additional_opts} \
            {params.model} \
            {input.file_dir} \
            > {output}.temp

        mv {output}.temp {output}
        """

# Rule dorado_demux demultiplexes a multiplexed run. The output path will
# contain multiple bam files in the format [kit]_barcode[barcode].bam (e.g.
# SQK-RPB004_barcode01.bam).
rule dorado_demux:
    input:
        EXP_DIR + "/{e}/dorado/calls.bam"
    output:
        directory = directory(EXP_DIR + "/{e}/dorado/demux"),
    threads:
        12
    resources:
        mem_mb = 64*1024,
        runtime = 2*24*60
    envmodules:
        "dorado/0.7.1"
    shell:
        """
        dorado demux \
            --no-classify \
            --no-trim \
            --threads {threads} \
            --output-dir {output.directory} \
            {input}
        """


# Rule dorado_demux_isolate_selected_bam moves the requested bam corresponding
# to a specific barcode to a different directory that can be accessed by
# subsequent rules. This rule essentially acts as an interface to the demux
# rule to describe that the output file names have a {kit} and {barcode}
# attribute.
rule dorado_demux_isolate_selected_bam:
    input:
        EXP_DIR + "/{e}/dorado/demux/",
    output:
        EXP_DIR + "/{e}/dorado/demux_selected/{kit}_barcode{b}.bam",
    shell:
        """
        mv {input}/{wildcards.kit}_barcode{wildcards.b}.bam {output}
        """


# dorado_bam_from_basecalling identifies and returns the path to the
# basecalled data for a sample. For barcoded samples the path contains two
# extra levels corresponding to the barcode.
def dorado_bam_from_basecalling(wilds):
    s = samples[wilds.s]

    if s.is_barcoded():
        # e.g. SQK-RPB004_barcode01.bam
        return os.path.join(EXP_DIR, s.parent_exp, "dorado",
                            "demux_selected", s.kit + '_barcode' + s.barcode + '.bam')

    return os.path.join(EXP_DIR, s.parent_exp, "dorado",
                        "calls.bam")


# Rule dorado_get_fastq_from_basecalled_bam_for_sample finds and converts the
# basecalled bam file corresponding to the requested sample {s} to fastq.
rule dorado_get_fastq_from_basecalled_bam_for_sample:
    input: 
        dorado_bam_from_basecalling
    output: 
        RAW_DIR + "/fastq/{s}_reads.fastq.gz"
    threads: 
        10
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools fastq --threads {threads} {input} \
                | pigz \
                > {output}
        """


# dorado_pychopper_trim_orient_reads uses pychopper to identify and trim the
# ONT barcodes. It also orients the reads 5' to 3'. This is only used for the
# cDNA protocol.
rule dorado_pychopper_trim_orient_reads:
    input:
        "{prefix}.fastq.gz"
    output:
        stats_output = "{prefix}.pychop.stats.tsv",
        report = "{prefix}.pychop.report.pdf",
        rescued = "{prefix}.pychop.rescued.fastq.gz",
        unclass = "{prefix}.pychop.unclass.fastq.gz",
        trimmed = "{prefix}.pychop.trimmed.fastq.gz",
    threads: 8
    resources:
        mem_mb = 20*1024,
        runtime = 1*24*60,
        disk_mb = 20*1024
    conda:
        "../envs/pychopper.yml"
    shell:
        """
        pychopper \
            -S {output.stats_output} \
            -r {output.report} \
            -k PCS111 \
            -t {threads} \
            -u >(gzip -c > {output.unclass}) \
            -w >(gzip -c > {output.rescued}) \
            {input} \
            - \
            | gzip -c > {output.trimmed}
        """


# dorado_pychopper_merge_trimmed_rescued merges the rescued and trimmed reads
# from pychopper into a single file.
rule dorado_pychopper_merge_trimmed_rescued:
    input:
        rescued = "{prefix}.pychop.rescued.fastq.gz",
        trimmed = "{prefix}.pychop.trimmed.fastq.gz",
    output:
        "{prefix}.pychopped.fastq.gz"
    threads: 2
    resources:
        mem_mb = 2*1024,
        runtime = 24*60
    shell:
        """
        cat {input.trimmed} {input.rescued} > {output}
        """


# dorado_pychopper_path_to_stranded_fastq identifies and returns the proper
# stranded fastq file depending on whether PCR-cDNA (pychopper had to run) or
# dRNA-seq was run.
def dorado_pychopper_path_to_stranded_fastq(sample):
    s = sample

    if s.is_unstranded():
        return os.path.join(
            SAMPLES_DIR, s.name, "fastq", "reads.pychopped.fastq.gz")

    return os.path.join(
        SAMPLES_DIR, s.name, "fastq", "reads.fastq.gz")


localrules: dorado_rename_final_stranded_fastq


# dorado_rename_final_stranded_fastq simply selects the pychopped or
# non-pychopped (if already stranded) fastq file and copies it to the
# destination.
rule dorado_rename_final_stranded_fastq:
    input:
        lambda ws: dorado_pychopper_path_to_stranded_fastq(samples[ws.sample])
    output:
        RAW_DIR + "/fastq/{s}_reads.fastq.gz"
    shell:
        """
        mv {input} {output}
        """
