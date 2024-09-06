# sort_bam sorts a bam file.
rule sort_bam:
    input:
        "{filename}.bam",
    output:
        "{filename}.sorted.bam",
    threads: 20
    resources:
        mem_mb=30*1024,
        runtime=3*60,
        disk_mb=100*1024,
    conda:
        "../envs/samtools.yml"
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
        "{filename}.bam",
    output:
        "{filename}.bam.bai",
    threads: 20
    resources:
        mem_mb=5*1024,
        runtime=3*60
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools index -@ {threads} {input}
        """


# bamfile_flagstats outputs alignment statistics for alignments.
rule bamfile_flagstats:
    input:
        "{filename}.bam",
    output:
        "{filename}.bam.flagstats.txt",
    threads: 4
    resources:
        mem_mb=5*1024,
        runtime=3*60
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools flagstat -O tsv --threads {threads} {input} > {output}
        """
