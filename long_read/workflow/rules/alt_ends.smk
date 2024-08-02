# Define parameters. Will use defaults if cannot get them from config
ALT_ENDS_DIR = config.get('ALT_ENDS_DIR', 'alt_ends')
ALT_ENDS_RES = config.get('ALT_ENDS_RES', ALT_ENDS_DIR + '/results')


localrules: count_reads_close_to_alternative_poly_sites


rule alt_ends_all:
    input:
        expand(
                ALT_ENDS_RES + "/{s}/alt_polya_site_counts.bed",
                s = samples.keys()),


# count_reads_close_to_alternative_poly_sites counts number of reads that are
# within a given region surrounding alternative poly(A) sites given in BED.
rule count_reads_close_to_alternative_poly_sites:
    input:
        bed = DATA_DIR + "/ebola/genes.alternative_polya_sites.bed",
        bam = SAMPLES_DIR + "/{s}/align/reads.toEbolaGenome.sorted.bam",
    output:
        ALT_ENDS_RES + "/{s}/alt_polya_site_counts.bed"
    shell:
        """
        python3 workflow/scripts/alt_ends/count_ends_in_bed.py \
                --bed {input.bed} \
                --bam {input.bam} \
                --extend-start 5 \
                --extend-end 5 \
                > {output}
        """
