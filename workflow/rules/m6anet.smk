# Define parameters
M6ANET_DIR = config.get('M6ANET_DIR', 'm6anet')
M6ANET_RES = config.get('M6ANET_RES', M6ANET_DIR + '/results')


localrules: m6anet_all


rule m6anet_all	:
    input:
        expand(
               M6ANET_RES + "/{s}/toTranscriptome/data.site_proba.csv",
               s = [s.name for s in samples.values() if s.is_direct_rna()]),
        expand(
               M6ANET_RES + "/{s}/toEbolaGenome/data.site_proba.csv",
               s = [s.name for s in samples.values() if s.is_direct_rna() and "24h" in s.name])


rule m6anet_dataprep:
    input:
        ea = NANOPOLISH_RES + "/{s}/{target}.eventalign.txt"
    output:
        M6ANET_RES + "/{s}/{target}/eventalign.index"
    params:
        odir = lambda wilds, output: os.path.dirname(output[0]),
    threads: 10
    resources:
        runtime=1*24*60,
        mem_mb=50*1024,
    conda:
        "../envs/m6anet.yml"
    shell:
        """
        m6anet dataprep \
            --eventalign {input.ea} \
            --out_dir {params.odir} \
            --readcount_max 1000000 \
            --n_processes {threads}
        """


rule m6anet_runinference:
    input:
        M6ANET_RES + "{s}/{target}/eventalign.index"
    output:
        M6ANET_RES + "{s}/{target}/data.site_proba.csv",
        M6ANET_RES + "{s}/{target}/data.indiv_proba.csv",
    params:
        work_dir = lambda wilds, output: os.path.dirname(output[0]),
    threads: 8
    resources:
        runtime=4*24*60,
        mem_mb=50*1024,
    conda:
        "../envs/m6anet.yml"
    shell:
        """
        m6anet inference \
            --input_dir {params.work_dir} \
            --out_dir {params.work_dir} \
            --n_processes {threads}
        """
