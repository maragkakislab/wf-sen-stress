# Define parameters. Will use defaults if cannot get them from config
RELPOS_DIR = config.get('RELPOS_DIR', 'relpos')
RELPOS_RES = config.get('RELPOS_RES', RELPOS_DIR + '/results')


rule relpos_all:
    input:
        expand(
            RELPOS_RES + "/{s}/{target}_{rpos}-vs-{fpos}.{region}.pdf",
            s = samples.keys(), target = ['toEbolaGenome'],
            rpos = ['5p', '3p'], fpos = ['5p', '3p'],
            region = ['mRNA', 'CDS']),


rule relpos_ebola:
    input:
        bam = SAMPLES_DIR + "/{s}/align/reads.{target}.sorted.bam",
        bed = DATA_DIR + '/' + 'ebola' + "/genes.{region}.bed",
    output:
        plot = RELPOS_RES +
               "/{s}/{target}_{rpos}-vs-{fpos}.{region}.pdf",
        text = RELPOS_RES +
               "/{s}/{target}_{rpos}-vs-{fpos}.{region}.txt",
    conda:
        "../envs/biometatools.yml"
    resources:
        mem_mb=2*1024,
        runtime=2*60
    shell:
        """
        rel-pos-bam-vs-bed.py \
            -m {input.bam} \
            -b {input.bed} \
            --fpos {wildcards.fpos} \
            --rpos {wildcards.rpos} \
            --plot-feats \
            -o {output.plot} \
            > {output.text}
        """
