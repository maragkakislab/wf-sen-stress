
# guppy_barcode_options returns the additional options that are used in the
# guppy_basecaller to enable detection of barcodes.
def guppy_barcode_options(experiment):
    if experiment.is_barcoded:
        return '--detect_barcodes --barcode_kits ' + experiment.kit
    return ''

# Rule guppy_basecall runs the guppy basecaller.
# Finds the folder gerneated from nanopore run and retrieves the fast5 files used for basecalling.
# {e} folder is generated manually.
rule guppy_basecall:
    input:
        RAW_DIR + "/{e}/origin.txt"
    output:
        EXP_DIR + "/{e}/guppy/sequencing_summary.txt"
    log:
        EXP_DIR + "/{e}/log/guppy.log"
    params:
        fc = lambda wilds: experiments[wilds.e].fc,
        kit = lambda wilds: experiments[wilds.e].kit,
        barcode_opts = lambda wilds:
            guppy_barcode_options(experiments[wilds.e]),
    threads:
        8
    resources:
        gpu = 2,
        gpu_model = "[gpuv100x|gpua100]",
        mem_mb = 64*1024,
        runtime = 10*24*60
    envmodules:
        "guppy/6.1.2"
    shell:
        """
        guppy_basecaller \
            -x cuda:all \
            --flowcell {params.fc} \
            --kit {params.kit} \
            --records_per_fastq 0 \
            --u_substitution off \
            --trim_strategy none \
            --input_path {RAW_DIR}/{wildcards.e}/runs/ \
            --save_path {EXP_DIR}/{wildcards.e}/guppy/ \
            --recursive \
            --gpu_runners_per_device 1 \
            --num_callers {threads} \
            --chunks_per_runner 512 \
            --compress_fastq \
            --max_queued_reads 20000 \
            {params.barcode_opts} \
            &> {log}
        """


# Rule merge_fastqs merges output fastqs per barcode from the basecalling step
# into a single file.
rule merge_fastqs_per_barcode:
    input:
        EXP_DIR + "/{e}/guppy/sequencing_summary.txt",
    output:
        fastq = EXP_DIR + "/{e}/guppy/barcode/{b}/reads.fastq.gz",
    shell:
        """
        find \
            {EXP_DIR}/{wildcards.e}/guppy/pass/barcode{wildcards.b}/ \
            {EXP_DIR}/{wildcards.e}/guppy/fail/barcode{wildcards.b}/ \
            -name "fastq_runid_*.fastq.gz" \
            -print0 \
        | while read -d $'\\0' FILE; do \
            cat $FILE >> {output.fastq};\
            done
        """

# input_for_get_fastq_from_basecalling identifies and returns the directory
# path that contains the basecalled data for a sample. For barcoded samples
# the path contains two extra levels corresponding to the barcode.
def input_for_get_fastq_from_basecalling(wilds):
    s = samples[wilds.s]

    if s.is_barcoded():
        return os.path.join(
            EXP_DIR, s.parent_exp, "guppy", "barcode", s.barcode, "reads.fastq.gz")

    return os.path.join(
        EXP_DIR, s.parent_exp, "guppy", "reads.fastq.gz")


# Rule merge_logs merges all logs from the basecalling step.
rule merge_logs:
    input:
        summ = EXP_DIR + "/{expt}/guppy/sequencing_summary.txt",
    output:
        glog = EXP_DIR + "/{expt}/guppy/guppy_basecaller.log.gz"
    shell:
        """
        find {EXP_DIR}/{wildcards.expt}/guppy/ \
            -name "guppy_basecaller_log-*.log" \
            -print0 \
        | while read -d $'\\0' FILE; do \
            cat $FILE | pigz -c >> {output.glog};\
            done
        """


# Rule clean_guppy_logs cleans all guppy output logs.
rule clean_guppy_logs:
    input:
        merged_log = EXP_DIR + "/{expt}/guppy/guppy_basecaller.log.gz",
    output:
        EXP_DIR + "/{expt}/guppy/clean_guppy_logs_done"
    shell:
        """
        find {EXP_DIR}/{wildcards.expt}/guppy/ \
            -name "guppy_basecaller_log-*.log" \
            -delete
        touch {output}
        """


# Rule clean_guppy cleans all output fastqs and logs from the basecalling step.
rule clean_guppy_fastqs:
    input:
        merge_fastqs_done = EXP_DIR + "/{expt}/guppy/reads.fastq.gz",
    output:
        EXP_DIR + "/{expt}/guppy/clean_guppy_fastqs_done"
    shell:
        """
        find {EXP_DIR}/{wildcards.expt}/guppy/ \
            -name "fastq_runid_*.fastq.gz" \
            -delete
        touch {output}
        """


# Rule clean_guppy cleans all output fastqs and logs from the basecalling step.
rule clean_guppy_fastqs_for_barcode:
    input:
        merge_fastqs_done = EXP_DIR + "/{e}/guppy/barcode/{b}/reads.fastq.gz",
    output:
        EXP_DIR + "/{e}/guppy/clean_guppy_fastqs_barcode_{b}_done"
    shell:
        """
        find \
            {EXP_DIR}/{wildcards.e}/guppy/pass/barcode{wildcards.b}/ \
            {EXP_DIR}/{wildcards.e}/guppy/fail/barcode{wildcards.b}/ \
            -name "fastq_runid_*.fastq.gz" \
            -delete
        touch {output}
        """




# Rule get_fastq_from_basecalling creates a symbolic link to the basecalled
# fastq file.
rule get_fastq_from_basecalling:
    input: input_for_get_fastq_from_basecalling
    output: SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz"
    threads: 1
    shell:
        """
        mv {input} {output}
        """


# Rule sanitize_headers removes extraneous information from the fastq headers
# to save space in subsequent steps.
rule sanitize_headers:
    input:
        SAMPLES_DIR + "/{s}/fastq/reads.fastq.gz"

    output:
        SAMPLES_DIR + "/{s}/fastq/reads.sanitize.fastq.gz"
    threads: 12
    resources:
        mem_mb = 2*1024,
        runtime = 6*60
    shell:
        """
        zcat {input} \
            | {SCRIPTS_DIR}/fastq_sanitize_header.py --fastq - \
            | pigz \
            > {output}
        """


# pychopper_trim_orient_reads uses pychopper to identify and trim the ONT
# barcodes. It also orients the reads 5' to 3'. This is only used for the cDNA
# protocol.
rule pychopper_trim_orient_reads:
    input:
        SAMPLES_DIR + "/{s}/fastq/{prefix}.fastq.gz"
    output:
        report = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.report.pdf",
        rescued = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.rescued.fastq.gz",
        unclass = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.unclass.fastq.gz",
        trimmed = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.trimmed.fastq.gz",
    threads: 8
    resources:
        mem_mb = 50*1024,
        runtime = 1*24*60,
        disk_mb = 200*1024
    envmodules:
        "pychopper/2.7.1"
    shell:
        """
        pigz -p {threads} -d -c {input} \
            > $TMPDIR/input.fastq

        pychopper \
            -r {output.report} \
            -k PCS111 \
            -t {threads} \
            -u >(gzip -c > {output.unclass}) \
            -w >(gzip -c > {output.rescued}) \
            $TMPDIR/input.fastq \
            - \
            | gzip -c > {output.trimmed}
        """


# pychopper_trim_orient_reads uses pychopper to identify and trim the ONT
# barcodes. It also orients the reads 5' to 3'. This is only used for the cDNA
# protocol.
rule pychopper_merge_trimmed_rescued:
    input:
        rescued = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.rescued.fastq.gz",
        trimmed = SAMPLES_DIR + "/{s}/fastq/{prefix}.pychop.trimmed.fastq.gz",
    output:
        SAMPLES_DIR + "/{s}/fastq/{prefix}.pychopped.fastq.gz",
    threads: 2
    resources:
        mem_mb = 2*1024,
        runtime = 24*60
    shell:
        """
        cat {input.trimmed} {input.rescued} > {output}
        """
