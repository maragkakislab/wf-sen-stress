rule flair_bam2bed:
    input:
        SAMPLES_DIR + "/{sample}/align/reads.toGenome.sorted.bam"
    output:
        temp(FLAIR_RES + "/{sample}/reads.toGenome.sorted.bed")
    resources:
        mem_mb = 30*1024
    conda:
        "../envs/flair.yml"
    shell:
        """
        bam2Bed12 -i {input} > {output}
        """


rule flair_correct:
    input:
        bed = FLAIR_RES + "/{sample}/reads.toGenome.sorted.bed",
        genome = GENOME_FILE,
        gtf = GTF_FILE,
    output:
        FLAIR_RES + "/{sample}/reads_all_corrected.bed",
        FLAIR_RES + "/{sample}/reads_all_inconsistent.bed"
    params:
        out_prefix = FLAIR_RES + "/{sample}/reads",
    resources:
        mem_mb = 200*1024,
        runtime = 24*60,
        disk_mb = 30
    threads: 16
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair correct \
            --threads {threads} \
            --query {input.bed} \
            --gtf {input.gtf} \
            --genome {input.genome} \
            --output {params.out_prefix}
        """


rule flair_concatenate_bed_files:
    input:
        bed = expand(FLAIR_RES + "/{s}/reads_all_corrected.bed", s = samples.keys())
    output:
        temp(FLAIR_RES + "/all/reads_all_corrected.bed")
    resources:
        mem_mb = 20*1024,
        runtime = 12*60,
        disk_mb = 30
    threads: 2
    envmodules:
       "bedops/2.4.41"
    shell:
        """
        bedops -u {input.bed} > {output}
        """


def input_reads_for_flair_collapse(samples):
    files = []
    for s in samples.values():
        if s.is_unstranded():
            files.append(os.path.join(
                SAMPLES_DIR, s.name, "fastq", "reads.pychopped.fastq.gz"))

        files.append(os.path.join(
            SAMPLES_DIR, s.name, "fastq", "reads.fastq.gz"))
    return files


rule flair_collapse:
    input:
        bed = FLAIR_RES + "/all/reads_all_corrected.bed",
        reads = lambda ws: input_reads_for_flair_collapse(samples),
        genome = GENOME_FILE,
        gtf = GTF_FILE,
    output:
        bed = FLAIR_RES + "/all/reads.isoforms.bed",
        fa = FLAIR_RES + "/all/reads.isoforms.fa",
        gtf =  FLAIR_RES + "/all/reads.isoforms.gtf",
    params:
        out_prefix = FLAIR_RES + "/all/reads",
    resources:
        mem_mb = 256*1024,
        runtime = 2*24*60,
        disk_mb = 120
    threads: 40
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair collapse \
            --threads {threads} \
            --query {input.bed} \
            --genome {input.genome} \
            --reads {input.reads} \
            --gtf {input.gtf} \
            --output {params.out_prefix}
        """


rule flair_make_metadata:
    output:
        FLAIR_RES + "/{sgroup}/metadata.tsv"
    threads: 1
    run:
        with open(output[0], 'w') as out:
            meta_list = FLAIR_METADATA[wildcards.sgroup]
            for l in meta_list:
                out.write("\t".join(l) + '\n')


rule flair_quantify:
    input:
        fa =  FLAIR_RES + "/all/reads.isoforms.fa",
        meta =  FLAIR_RES + "/all/metadata.tsv"
    output:
        tsv = FLAIR_RES + "/all/reads.flair.quantify",
        temp_dir = temp(directory(FLAIR_RES + "/all/reads.flair.quantify.temp"))
    params:
        out_prefix = FLAIR_RES + "/all/reads.flair.quantify",
    resources:
        mem_mb = 120*1024,
        runtime = 12*60,
        disk_mb = 60
    threads: 40
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair quantify \
            --reads_manifest {input.meta} \
            --isoforms {input.fa} \
            --threads {threads} \
            --temp_dir {output.tsv} \
            --output {params.out_prefix}
        """


rule flair_diffexp:
    input:
        tsv = FLAIR_RES + "/all/reads.flair.quantify"
    output:
        touch(FLAIR_RES + "/all/reads.flair.diffexp/done")
    params:
        outdir = lambda wilds, output: os.path.dirname(output[0]),
    resources:
        mem_mb = 50*1024,
        runtime = 5*60,
        disk_mb = 60
    threads: 20
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair diffexp \
            --counts_matrix {input.tsv} \
            --threads {threads} \
            --out_dir {params.outdir} \
            --out_dir_force
        """


rule flair_diffsplice:
    input:
        tsv = FLAIR_RES + "/all/reads.flair.quantify",
        bed = FLAIR_RES + "/all/reads.isoforms.bed"
    output:
        FLAIR_RES + "/all/reads.flair.diffsplice/diffsplice.alt3.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice/diffsplice.alt5.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice/diffsplice.es.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice/diffsplice.ir.events.quant.tsv",
    params:
        outdir = lambda wilds, output: os.path.dirname(output[0]),
    resources:
        mem_mb = 50*1024,
        runtime = 5*60,
        disk_mb = 60
    threads: 20
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair diffsplice \
            --isoforms {input.bed} \
            --counts_matrix {input.tsv} \
            --threads {threads} \
            --test \
            --out_dir {params.outdir} \
            --out_dir_force
        """
