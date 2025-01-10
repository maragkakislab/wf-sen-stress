# align_reads_to_genome aligns the input reads to the genome.
ASSEMBLY = config["ASSEMBLY"]

rule align_reads_to_genome:
    input:
        fastq = RAW_DIR + "/fastq/{e}/{s}_reads.fastq.gz",
        genome = DATA_DIR + /ASSEMBLY + "/genome/genome.fa",
    output:
        temp(SAMPLES_DIR + "/{s}/reads.toGenome.bam"),
    threads: 50
    resources:
        mem_mb = 100*1024,
        runtime = 2*24*60
    conda:
        "../envs/minimap2.yml"
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
        fastq = RAW_DIR + "/fastq/{e}/{s}_reads.fastq.gz",
        transcriptome = DATA_DIR + "/hg38/transcriptome/transcripts.fa",
    output:
        temp(SAMPLES_DIR + "/{sample}/reads.toTranscriptome.bam"),
    threads: 50
    resources:
        mem_mb = 100*1024,
        runtime = 2*24*60
    conda:
        "../envs/minimap2.yml"
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
        SAMPLES_DIR + "/{e}/{s}/{prefix}.bam",
    output:
        SAMPLES_DIR + "/{e}/{s}/{prefix}.sorted.bam",
    threads: 20
    resources:
        mem_mb=30*1024,
        runtime=2*60,
        disk_mb=100*1024,
    conda:
        "../envs/minimap2.yml"
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
        SAMPLES_DIR + "/{e}/{s}/{prefix}.bam",
    output:
        SAMPLES_DIR + "/{e}/{s}/{prefix}.bam.bai",
    threads: 40
    resources:
        mem_mb=5*1024,
        runtime=5*60
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        samtools index -@ {threads} {input}
        """



# bamfile_flagstats outputs alignment statistics for alignments.
rule bamfile_flagstats:
    input:
        SAMPLES_DIR + "/{e}/{s}/{prefix}.bam",
    output:
        SAMPLES_DIR + "/{e}/{s}/{prefix}.bam.flagstats.txt",
    threads: 4
    resources:
        mem_mb=5*1024,
        runtime=3*60
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        samtools flagstat -O tsv --threads {threads} {input} > {output}
        """

# count_aligned_reads_per_transcript counts the reads aligned on each
# transcript.
rule count_aligned_reads_per_transcript:
    input:
        aligned = SAMPLES_DIR + "{e}/{s}/reads.toTranscriptome.bam",
        transcript_tab = DATA_DIR + "/hg38/transcript-gene.tab",
    output:
        ANALYSIS_DIR + "/counts/{s}/reads.toTranscriptome.txt",
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        {SCRIPTS_DIR}/sam_per_ref_count_statistics.py \
            --ifile {input.aligned} \
            --ref-col-name transcript \
            --cnt-col-name count \
            --opt-col-name sample \
            --opt-col-val {wildcards.sample} \
            | table-join.py \
                --table1 - \
                --table2 {input.transcript_tab} \
                --key1 transcript \
                --key2 transcript \
                > {output}
        """

####Aggregate transcriptome aligned reads to gene level

rule transcriptome_to_gene_level:
    input:
        ANALYSIS_DIR + "/counts/{s}/reads.toTranscriptome.txt",
    output:
        ANALYSIS_DIR + "/counts/{s}/reads.toGenome.txt",
    run:
        import pandas as pd
        f = input[0]
        x = pd.read_csv(f, sep = '\t')
        y = x.loc[:,["sample","count","gene"]].groupby(["gene","sample"],as_index = False).sum()
        y.to_csv(output[0], sep='\t')