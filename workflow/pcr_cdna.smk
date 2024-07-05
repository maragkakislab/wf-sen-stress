import sys
import glob


class sample:
    def __init__(self, name, kit, fc, parent_exp, barcode=None):
        self.name = name
        self.kit = kit
        self.fc = fc
        self.parent_exp = parent_exp
        self.barcode = barcode

    def is_barcoded(self):
        if self.barcode is None:
            return False
        return True

    def is_direct_rna(self):
        if self.kit == 'SQK-RNA002':
            return True
        return False

    def is_unstranded(self):
        if self.kit in ['SQK-PCB111-24', 'SQK-PCS111']:
            return True
        return False


class experiment:
    def __init__(self, name, kit, fc, is_barcoded):
        self.name = name
        self.kit = kit
        self.fc = fc
        self.is_barcoded = is_barcoded


# Load config
if "--configfile" in sys.argv:
    i = sys.argv.index("--configfile")
    config_path = sys.argv[i + 1]
    configfile: config_path
else:
    config_path = "config/config.yml"
    configfile: config_path


# Set variables
# EXP_DIR must be generated manually and contain experiment files from the promethion
ANALYSIS_DIR = config['ANALYSIS_DIR']
DATA_DIR = config['DATA_DIR']
RAW_DIR = config['RAW_DIR']
EXP_DIR = config['EXP_DIR']
IGV_DIR = config['IGV_DIR']
SAMPLES_DIR = config['SAMPLES_DIR']
SCRIPTS_DIR = config['SCRIPTS_DIR']
TRASH_DIR = config['TRASH_DIR']

RSYNC_PATH = config['RSYNC_PATH']



# Create a dictionary with samples and experiments.
# Takes a list of samples provided in the Config file under SAMPLE DATA and turns them into the sample class.
# These sample names will be wildcard {s}
# Generates a list of experiments in parrallel
samples = {}
experiments = {}
for d in config['SAMPLE_DATA']:
    s = sample(*d)
    samples[s.name] = s
    experiments[s.parent_exp] = experiment(s.parent_exp, s.kit, s.fc,
                                           s.is_barcoded())


##### include rules #####
include: "rules/get_data.smk"
include: "rules/annotation.smk"

# Define rules that require minimal resources and can run on the same node
# where the actual snakemake is running
localrules: run_all, merge_fastqs, merge_fastqs_per_barcode, merge_logs,
            clean_guppy_logs, get_fastq_from_basecalling,
            symlink_alignments_for_igv, clean_guppy_fastqs,
            clean_guppy_fastqs_for_barcode, aggregate

# Rule run_all collects all outputs to force execution of the whole pipeline.
# Identified files will be produced
rule run_all:
    input:
        # Alignments BAM/BAI files
        expand(
            SAMPLES_DIR + "/{s}/align/reads.{target}.sorted.{sufx}",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome'],
            sufx=['bam', 'bam.bai']),

        # Flagstats
        expand(
            SAMPLES_DIR + "/{s}/align/reads.{target}.sorted.bam.flagstats.txt",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome']),

        # Counts per transcript/gene
        expand(
            ANALYSIS_DIR +
            "/counts/{s}/reads.toTranscriptome.txt",
            s=samples.keys(),
            ),

        
         # Counts per gene
        expand(
            ANALYSIS_DIR +
            "/counts/{s}/reads.toGenome.txt",
            s=samples.keys(),
            ),
 
        # collated tables
        ANALYSIS_DIR + "/all_genome_counts.txt",
        ANALYSIS_DIR + "/all_transcriptome_counts.txt",


        # # DESpline
        # ANALYSIS_DIR + "/spline/DE_Spline_res.rds",



        # # Relative positions
        # rules.relpos_all.input,

        # # Files for IGV
        # expand(
        #     IGV_DIR + "/align/{s}-reads.{target}.sorted.{sufx}",
        #     s=samples.keys(),
        #     target=['toGenome', 'toTranscriptome', 'toEbolaGenome'],
        #     sufx=['bam', 'bam.bai']),

        # # Run fastqc for quality control
        # expand(
        #     SAMPLES_DIR + "/{s}/qc/reads_fastqc.html",
        #     s=samples.keys()),

        # Cleanup guppy files
        expand(
            EXP_DIR + "/{e}/guppy/clean_guppy_logs_done",
            e=experiments.keys()),

        # Cleanup guppy fastq files
        expand(
            EXP_DIR + "/{e}/guppy/clean_guppy_fastqs_done",
            e=[e.name for e in experiments.values() if not e.is_barcoded]),
        expand(
            EXP_DIR + "/{s.parent_exp}/guppy/clean_guppy_fastqs_barcode_{s.barcode}_done",
            s=[s for s in samples.values() if s.barcode is not None]),



############################################### Rules


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


# input_for_alignment rules identifies and returns the input fastq file for
# the alignment step. For unstranded samples the file will pass through
# pychopper to orient and remove adapters.
def input_for_alignment_rules(sample, fastq_prefix):
    s = sample

    if s.is_unstranded():
        return os.path.join(
            SAMPLES_DIR, s.name, "fastq", fastq_prefix + ".pychopped.fastq.gz")

    return os.path.join(
        SAMPLES_DIR, s.name, "fastq", fastq_prefix + ".fastq.gz")


# align_reads_to_genome aligns the input reads to the genome.
rule align_reads_to_genome:
    input:
        fastq = lambda ws: input_for_alignment_rules(samples[ws.sample], ws.prefix),
        genome = DATA_DIR + "/hg38/genome/genome.fa",
    output:
        SAMPLES_DIR + "/{sample}/align/{prefix}.toGenome.bam",
    threads: 50
    resources:
        mem_mb = 100*1024,
        runtime = 2*24*60
    conda:
        "envs/minimap2.yml"
    envmodules:
        "samtools/1.15.1",
        "minimap2/2.24"
    shell:
        """
        minimap2 \
                -a \
                -x splice \
                -k 12 \
                -u b \
                --MD \
                --sam-hit-only \
                -t {threads} \
                --secondary=no \
                {input.genome}\
                {input.fastq} \
                    | grep -v "SA:Z:" \
                    | samtools view -b -F 256 - \
                    > {output}
        """


# align_reads_to_transcriptome: aligns the input reads to the transcriptome.
rule align_reads_to_transcriptome:
    input:
        fastq =  lambda ws: input_for_alignment_rules(samples[ws.sample], ws.prefix),
        transcriptome = DATA_DIR + "/hg38/t7_transcripts.fa",
    output:
        SAMPLES_DIR + "/{sample}/align/{prefix}.toTranscriptome.bam",
    threads: 50
    resources:
        mem_mb = 100*1024,
        runtime = 2*24*60
    conda:
        "envs/minimap2.yml"
    envmodules:
        "samtools/1.15.1",
        "minimap2/2.24"
    shell:
        """
        minimap2 \
                -a \
                -x map-ont \
                -k 12 \
                -u f \
                -t {threads} \
                --secondary=no \
                {input.transcriptome}\
                {input.fastq} \
                    | grep -v "SA:Z:" \
                    | samtools view -b -F 256 - \
                    > {output}
        """

# sort_bam sorts a bam file.
rule sort_bam:
    input:
        SAMPLES_DIR + "/{sample}/align/{filename}.bam",
    output:
        SAMPLES_DIR + "/{sample}/align/{filename}.sorted.bam",
    threads: 20
    resources:
        mem_mb=30*1024,
        runtime=2*60,
        disk_mb=100*1024,
    conda:
        "envs/minimap2.yml"
    envmodules:
        "samtools/1.15.1",
    shell:
        """
            samtools sort \
                --threads {threads} \
                -T /lscratch/$SLURM_JOB_ID \
                -o {output} \
                {input}
        """

# index_bam indexes a bam file
rule index_bam:
    input:
        SAMPLES_DIR + "/{sample}/align/{filename}.bam",
    output:
        SAMPLES_DIR + "/{sample}/align/{filename}.bam.bai",
    threads: 40
    resources:
        mem_mb=5*1024,
        runtime=5*60
    envmodules:
        "samtools/1.15.1"
    shell:
        """
        samtools index -@ {threads} {input}
        """

# symlink_alignments_for_igv creates symbolic links for bam files to be used
# with IGV. The links contain the sample name in the filename to simplify
# loading in IGV.
rule symlink_alignments_for_igv:
    input:
        SAMPLES_DIR + "/{sample}/align/{filename}",
    output:
        "igv/align/{sample}-{filename}",
    threads: 1
    shell:
        """
        ln -sf ../../{input} {output}
        """


# bamfile_flagstats outputs alignment statistics for alignments.
rule bamfile_flagstats:
    input:
        SAMPLES_DIR + "/{sample}/align/{filename}.bam",
    output:
        SAMPLES_DIR + "/{sample}/align/{filename}.bam.flagstats.txt",
    threads: 4
    resources:
        mem_mb=5*1024,
        runtime=3*60
    envmodules:
        "samtools/1.15.1"
    shell:
        """
        samtools flagstat -O tsv --threads {threads} {input} > {output}
        """

# count_aligned_reads_per_transcript counts the reads aligned on each
# transcript.
rule count_aligned_reads_per_transcript:
    input:
        aligned = SAMPLES_DIR + "/{s}/align/{prefix}.toTranscriptome.sorted.bam",
        transcript_tab = DATA_DIR + "/hg38/transcript-gene.tab",
    output:
        ANALYSIS_DIR + "/counts/{s}/{prefix}.toTranscriptome.txt",
    conda:
        "envs/bam.yml" 
    shell:
        """
        {SCRIPTS_DIR}/sam_per_ref_count_statistics.py \
            --ifile {input.aligned} \
            --ref-col-name transcript \
            --cnt-col-name count \
            --opt-col-name sample \
            --opt-col-val {wildcards.s} \
            | table-join.py \
                --table1 - \
                --table2 {input.transcript_tab} \
                --key1 transcript \
                --key2 transcript \
                > {output}
        """


####Aggregate transcriptome aligned reads to genome aligned reads

rule aggregate:
    input:
        ANALYSIS_DIR + "/counts/{s}/{prefix}.toTranscriptome.txt",
    output:
        ANALYSIS_DIR + "/counts/{s}/{prefix}.toGenome.txt",
    run:
        import pandas as pd
        f = input[0]
        x = pd.read_csv(f, sep = '\t')
        y = x.loc[:,["sample","count","gene"]].groupby(["gene","sample"],as_index = False).sum()
        y.to_csv(output[0], sep='\t')


# ###### Generate a single dataframe

rule input_for_DE_spline:
    input:
         expand(ANALYSIS_DIR + "/counts/{s}/reads.toGenome.txt", s=samples.keys()),
    output:
        ANALYSIS_DIR + "/all_genome_counts.txt",
    run:
        import pandas as pd
        dfs = [pd.read_csv(file_path, delimiter='\t') for file_path in input]
        combined_df = pd.concat(dfs, axis=0, ignore_index=True)
        wide_df = combined_df.pivot(index= "gene", columns= 'sample', values = 'count')
        wide_df.to_csv(output[0], sep = '\t')


rule transcriptome_combined:
    input:
         expand(ANALYSIS_DIR + "/counts/{s}/reads.toTranscriptome.txt", s=samples.keys()),
    output:
        ANALYSIS_DIR + "/all_transcriptome_counts.txt",
    run:
        import pandas as pd
        dfs = [pd.read_csv(file_path, delimiter='\t') for file_path in input]
        combined_df = pd.concat(dfs, axis=0, ignore_index=True)
        wide_df = combined_df.pivot(index= "transcript", columns= 'sample', values = 'count')
        wide_df.to_csv(output[0], sep = '\t')





# # # DESpline_R analysis
# rule DE_spline:
#     input:
#         counts =  "analysis/all_genome_counts.txt",
#         metadata = "data/metadata.txt",
#     output:
#         ANALYSIS_DIR + "/spline/DE_Spline_res.rds",
#     envmodules:
#         "R/4.3"
#     shell:
#         """
#         {SCRIPTS_DIR}/splinegroupR/src/DESeq_Spline.R \
#         -c {input.counts} \
#         -m {input.metadata} \
#         -o analysis/spline \
#         """




# # fastqc_fastq runs fastqc for a fastq file
# rule fastqc_fastq:
#     input:
#         SAMPLES_DIR + "/{s}/fastq/{prefix}.fastq.gz",
#     output:
#         SAMPLES_DIR + "/{s}/qc/{prefix}_fastqc.html",
#     params:
#         outdir = lambda wilds, output: os.path.dirname(output[0]),
#         dummy_threads = 20
#     threads: 2
#     resources:
#         mem_mb=int(20*250*1.5), # dummy_threads * 250M(see below) + 50% extra
#         runtime=24*60
#     envmodules:
#         "fastqc/0.11.9"
#     shell:
#         """
#             # Memory limit for FastQC is defined as threads*250M. Although
#             # it's not parallelized, below we use -t dummy_threads to
#             # indirectly increase the memory limit.
#             fastqc \
#                 -o {params.outdir} \
#                 -t {params.dummy_threads} \
#                 {input}
#         """





# ################################################################################
# # fastq_len_distro calculates the length distribution for a fastq file.
# ################################################################################
# rule fastq_len_distro:
#     input:
#         SAMPLES_DIR + "/{sample}/fastq/reads.sanitize.pychopper.fastq.gz",
#     output:
#         tab = SAMPLES_DIR + "/{sample}/qc/reads.sanitize.pychopper.fastq_lendistro.tab",
#         fig = SAMPLES_DIR + "/{sample}/qc/reads.sanitize.pychopper.fastq_lendistro.pdf",
#     threads: 40
#     resources: mem_mb=100*1024, runtime=24*60
#     shell:
#         """
#             /data/Maragkakislab/darsa/workspace/dev/go/bin/fastq-len-distro \
#                 --skip-zeros \
#                 {input} \
#                 | /data/Maragkakislab/darsa/workspace/dev/tabletools/table-paste-col/table-paste-col.py \
#                  --table - --col-name sample --col-val null \
#                 > {output.tab}

#             ./dev/plot-len-distro.R \
#                 --ifile {output.tab} \
#                 --figfile {output.fig}
#         """

# ################################################################################
# # count_fastq_reads: counts of the reads per sample for multiple files
# ################################################################################
# # NOTE: OPTIMIZE FOR single filename

# rule count_fastq_reads:
#     input:
#         fq_in_sanit = SAMPLES_DIR + "/{sample}/fastq/reads.sanitize.fastq.gz",
#         fq_in_pychop = SAMPLES_DIR + "/{sample}/fastq/reads.sanitize.pychopper.fastq.gz",
#         fq_in_rescued = SAMPLES_DIR + "/{sample}/fastq/rescued.fq.gz",
#         fq_in_unused = SAMPLES_DIR + "/{sample}/fastq/unclassified.fq.gz",
#         # SAMPLES_DIR + "/{sample}/fastq/{fastqfile}",

#     output:
#         fq_out_sanit = resdir_analysis + "/counts/results/{sample}/fastq/counts_reads.sanitize.fastq.gz.txt",
#         fq_out_pychop = resdir_analysis + "/counts/results/{sample}/fastq/counts_reads.sanitize.pychopper.fastq.txt",
#         fq_out_rescued = resdir_analysis + "/counts/results/{sample}/fastq/counts_rescued.fq.txt",
#         fq_out_unused = resdir_analysis + "/counts/results/{sample}/fastq/counts_unclassified.fq.txt",
#         # "resdir_counts" + "/{sample}/counts_{fastqfile}.txt"
#     threads: 4
#     resources: mem_mb=10*1024, runtime=2*60
#     shell:
#         """
#             awk '{{s++}}END{{print s/4}}' {input.fq_in_sanit} > {output.fq_out_sanit}
#             awk '{{s++}}END{{print s/4}}' {input.fq_in_pychop} > {output.fq_out_pychop}
#             awk '{{s++}}END{{print s/4}}' {input.fq_in_rescued} > {output.fq_out_rescued}
#             awk '{{s++}}END{{print s/4}}' {input.fq_in_unused} > {output.fq_out_unused}
#         """



# ################################################################################
# # count_aligned_reads: counts of the reads per sample for multiple files
# ################################################################################
# rule count_aligned_reads:
#     input:
#         SAMPLES_DIR + "/{sample}/align/reads.sanitize.pychopper.toTranscriptome.sorted.bam",
#     output:
#         resdir_analysis + "/counts/results/{sample}/align/read.counts.txt"
#     shell:
#         """
#         analysis/counts/src/sam_per_ref_count_statistics.py \
#             --ifile {input} \
#             --ref-col-name transcript \
#             --cnt-col-name count \
#             --opt-col-name sample \
#             --opt-col-val {wildcards.sample} \
#             | /data/Maragkakislab/darsa/workspace/dev/tabletools/table-join/table-join.py \
#                 --table1 - \
#                 --table2 /data/Maragkakislab/darsa/workspace/projects/indi/data/hg38/transcript-gene.tab \
#                 --key1 transcript \
#                 --key2 transcript \
#                 > {output}
#         """
# ################################################################################
# # quantile_norm_transcripts: Quantile normalisation of Transctipt level
# ################################################################################
# rule quantile_norm_transcripts:
#     input:
#         resdir_analysis + "/counts/results/{sample}/align/read.counts.txt",
#     output:
#         resdir_analysis + "/counts/results/{sample}/align/tid_qnorm.txt"
#     shell:
#         """
#         analysis/counts/src/quantile_normalization.py \
#             --ifile {input} \
#             --quantile 0.90 \
#             --norm-col count \
#             --new-col-name qnorm90 \
# 			| /data/Maragkakislab/darsa/workspace/dev/tabletools/table-cut-columns/table-cut-columns.py \
# 				--table - \
# 				--col-name transcript sample qnorm90 \
# 				> {output}
#         """


# ################################################################################
# # quantile_normalization_gene: Quantile normalisation  gene level
# ################################################################################
# rule quantile_norm_gene:
#     input:
#         resdir_analysis + "/counts/results/{sample}/align/read.counts.txt",
#     output:
#         resdir_analysis + "/counts/results/{sample}/align/gid_qnorm.txt"
#     shell:
#         """
#         /data/Maragkakislab/darsa/workspace/dev/tabletools/table-group-summarize/table-group-summarize.py \
# 			--table {input} \
# 			--groupby gene sample \
# 			--summarize count \
# 			--func sum \
#             | analysis/counts/src/quantile_normalization.py \
#                 --ifile - \
#                 --quantile 0.90 \
#                 --norm-col count_sum \
#                 --new-col-name qnorm90 \
#     			| /data/Maragkakislab/darsa/workspace/dev/tabletools/table-cut-columns/table-cut-columns.py \
#     				--table - \
#     				--col-name gene sample qnorm90 \
#     				> {output}
#         """

# ################################################################################
# # get_raw_counts: get raw cpounts per sample
# ################################################################################
# rule get_raw_counts:
#     input:
#         resdir_analysis + "/counts/results/{sample}/align/read.counts.txt",
#     output:
#         gid_counts = resdir_analysis + "/counts/results/{sample}/align/gid_raw_counts.txt",
#         tid_counts = resdir_analysis + "/counts/results/{sample}/align/tid_raw_counts.txt",
#     shell:
#         """
#         /data/Maragkakislab/darsa/workspace/dev/tabletools/table-group-summarize/table-group-summarize.py \
# 			--table {input} \
# 			--groupby gene sample \
# 			--summarize count \
# 			--func sum \
# 			| /data/Maragkakislab/darsa/workspace/dev/tabletools/table-cut-columns/table-cut-columns.py \
# 				--table - \
# 				--col-name gene count_sum sample \
# 				> {output.gid_counts}

# 		/data/Maragkakislab/darsa/workspace/dev/tabletools/table-cut-columns/table-cut-columns.py \
# 			--table {input} \
# 			--col-name transcript count sample \
# 			> {output.tid_counts}
#         """
# ################################################################################
# # combined_raw_counts: get raw cpounts per sample
# ################################################################################
# rule combined_raw_counts:
#     input:
#         in_gid_counts = expand(resdir_analysis + "/counts/results/{sample}/align/gid_raw_counts.txt", sample = samples),
#         in_tid_counts = expand(resdir_analysis + "/counts/results/{sample}/align/tid_raw_counts.txt", sample = samples),
#     output:
#         out_gid_counts = resdir_analysis + "/PCA/results/raw_data/all_gid_raw_counts.txt",
#         out_tid_counts = resdir_analysis + "/PCA/results/raw_data/all_tid_raw_counts.txt",
#     shell:
#         """
#         awk 'FNR>1 || NR==1' {input.in_gid_counts} > {output.out_gid_counts}
#         awk 'FNR>1 || NR==1' {input.in_tid_counts} > {output.out_tid_counts}
#         """


# ################################################################################
# # combined_raw_counts: get raw cpounts per sample
# ################################################################################
# rule raw_counts_longformat:
#     input:
#         in_gid_counts = resdir_analysis + "/PCA/results/raw_data/all_gid_raw_counts.txt",
#         in_tid_counts = resdir_analysis + "/PCA/results/raw_data/all_tid_raw_counts.txt",
#     output:
#         out_gid_colwise = resdir_analysis + "/PCA/results/raw_data/all_gid_colwise_raw_counts.txt",
#         out_tid_colwise = resdir_analysis + "/PCA/results/raw_data/all_tid_colwise_raw_counts.txt",
#     shell:
#         """
#         analysis/PCA/src/table_long_to_wide.py \
#         	--ifile {input.in_gid_counts} \
#         	--index gene \
#         	--columns sample \
#         	--values count_sum \
#         	> {output.out_gid_colwise}

#         analysis/PCA/src/table_long_to_wide.py \
#         	--ifile {input.in_tid_counts} \
#         	--index transcript \
#         	--columns sample \
#         	--values count \
#         	> {output.out_tid_colwise}
#         """

# ################################################################################
# ########################### E.O.F ##############################################
# # ###############################################################################
# # deseq2_normalisation:
# # ###############################################################################
# # rule deseq2_normalisation:
# #     input:
# #         in_gid_colwise = resdir_analysis + "/PCA/results/raw_data/all_gid_colwise_raw_counts.txt",
# #         # in_tid_counts = resdir_analysis + "/PCA/results/raw_data/all_tid_colwise_raw_counts.txt",
# #     output:
# #         out_gid_counts = resdir_analysis + "/PCA/results/deseq2_normalised__PCA.pdf",
# #         # out_tid_counts = resdir_analysis + "/PCA/results/raw_data/colwise_all_tid_raw_counts.txt",
# #     shell:
# #
# # shell:
# #     """
# #     module load R/3.6
# #     ./src/deseq2_normalization.R \
# #     	{input.in_gid_colwise} \
# #     	resdir_analysis + "/PCA/metadata.txt \
# #     	tab \
# #     	{output.in_gid_colwise}
# #     """

# ################################################################################
# ############################################################################

# # echo ">>> Convert long format to wide format <<<"
# # ############################################################################
# # ./src/table_long_to_wide.py \
# # 	--ifile "results/geneid/gid_combinedfile.tab" \
# # 	--index gene \
# # 	--columns sample \
# # 	--values count \
# # 	> "results/geneid/gid_rawcounts_col_wise.tab"
# #
# # ############################################################################
# # echo ">>> RUN DESeq2 analysis <<<"
# # ############################################################################
# # mkdir -p results/geneid_TRY/
# # ./src/deseq2_normalization.R \
# # 	results/geneid/gid_rawcounts_col_wise.tab \
# # 	results/metadata.txt \
# # 	tab \
# # 	results/geneid/deseq2_normalised_


# # softlink_fast5 creates symbolic link pointing to the directory with the raw
# # fast5 files
# rule softlink_fast5:
#     input:
#         origin=lambda wildcards: glob.glob(
#                 EXP_DIR + "/"+('[0-9]'*8)+"_{s}/origin.txt".format(s=wildcards.s)),
#     output:
#         SAMPLES_DIR + "/{s}/fast5linked"
#     params:
#         fast5dir=lambda wildcards: "../../" + glob.glob(
#                 EXP_DIR + "/" + ('[0-9]'*8) + "_{s}/".format(s=wildcards.s))[0] + "runs"
#     shell:
#         """
#         ln -s {params.fast5dir} {SAMPLES_DIR}/{wildcards.s}/fast5
#         touch {output}
#         """
