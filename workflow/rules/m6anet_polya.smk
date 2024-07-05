# Define parameters. Will use defaults if cannot get them from config
M6ANET_POLYA_DIR = config.get('M6ANET_POLYA_DIR', 'm6anet_polya')
M6ANET_POLYA_RES = config.get('M6ANET_POLYA_RES', M6ANET_POLYA_DIR + '/results')


localrules: m6anet_polya_all


rule m6anet_polya_all:
    input:
        expand(
                M6ANET_POLYA_RES + "/{s}/toEbolaGenome/polya_len_dist_per_mod_ratio.pdf",
                s = [s.name for s in samples.values() if s.is_direct_rna() and "24h" in s.name])


# plot_polya_len_dist_per_mod_ratio plots the polyA length distribution
# stratified by the modification ratio of the reads..
rule plot_polya_len_dist_per_mod_ratio:
    input:
        polya = NANOPOLISH_RES + "/{s}/ebola.polya.filtered.tab",
        bed = DATA_DIR + "/ebola/genes.mRNA.bed",
        bam = SAMPLES_DIR + "/{s}/align/reads.toEbolaGenome.sorted.bam",
        summary = NANOPOLISH_RES + "/{s}/toEbolaGenome.summary.txt",
        indiv_prob = M6ANET_RES + "/{s}/toEbolaGenome/data.indiv_proba.csv",
    output:
        M6ANET_POLYA_RES + "/{s}/toEbolaGenome/polya_len_dist_per_mod_ratio.pdf"
    conda:
        "../envs/python.yml"
    shell:
        """
        python3 workflow/scripts/m6anet_polya/polya_len_dist_per_mod_ratio.py \
                --summary {input.summary} \
                --indiv-prob {input.indiv_prob} \
                --bed {input.bed} \
                --bam {input.bam} \
                {input.polya} \
                {output}
        """
