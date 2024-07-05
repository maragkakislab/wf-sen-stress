import os


localrules: add_start_stop_codon_entries_to_gff, convert_gff_to_gtf,
            convert_gff_to_bed, filter_bed_for_feature


rule add_start_stop_codon_entries_to_gff:
    input:
        gff = "{filepath}/genes.gff3",
        genome = "{filepath}/genome.fa",
    output:
        "{filepath}/genes.corrected.gff3"
    conda:
        "../envs/agat.yml"
    shell:
        """
        agat_sp_add_start_and_stop.pl \
                --gff {input.gff} \
                --fasta {input.genome} \
                --out {output}
        """


rule convert_gff_to_gtf:
    input:
        "{filepath}/genes.corrected.gff3"
    output:
        "{filepath}/genes.corrected.gtf"
    conda:
        "../envs/agat.yml"
    shell:
        """
        agat_convert_sp_gff2gtf.pl --gff {input} -o {output}
        """


rule convert_gff_to_bed:
    input:
        "{filepath}/genes.gff3"
    output:
        "{filepath}/genes.bed"
    conda:
        "../envs/agat.yml"
    shell:
        """
        gff2bed < {input} > {output}
        """


rule filter_bed_for_feature:
    input:
        "{filepath}/genes.bed"
    output:
        "{filepath}/genes.{feature}.bed"
    shell:
        """
        grep -P "\t{wildcards.feature}\t" {input} > {output}
        """
