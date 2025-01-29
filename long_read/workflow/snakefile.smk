import sys
import glob
import pandas as pd
import re
# Load config
if "--configfile" in sys.argv:
    i = sys.argv.index("--configfile")
    config_path = sys.argv[i + 1]
    configfile: config_path
else:
    config_path ="config/config.yml"
    configfile: config_path

# Set variables
# EXP_DIR must be generated manually and contain experiment files from the promethion
ANALYSIS_DIR = config['ANALYSIS_DIR']
DATA_DIR = config['DATA_DIR']
RAW_DIR = config['RAW_DIR']
EXP_DIR = config['EXP_DIR']
SAMPLES_DIR = config['SAMPLES_DIR']
SCRIPTS_DIR = config['SCRIPTS_DIR']

# Create classes
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
        if self.kit in ['SQK-PCB111-24', 'SQK-PCS114-24']:
            return True
        return False


class experiment:
    def __init__(self, name, kit, fc, is_barcoded):
        self.name = name
        self.kit = kit
        self.fc = fc
        self.is_barcoded = is_barcoded


# Create a dictionary with samples and experiments.
# Takes a list of samples provided in the Config file under SAMPLE DATA and turns them into the sample class.
# These sample names will be wildcard {s} or {sample}
# Generates a list of experiments in parrallel
samples = {}
experiments = {}
for d in config['SAMPLE_DATA']:
    s = sample(*d)
    samples[s.name] = s
    experiments[s.parent_exp] = experiment(s.parent_exp, s.kit, s.fc,
                                           s.is_barcoded())
comparisons = config['COMPARISONS']



##### include rules #####
include: "rules/dorado.smk"
include: "rules/minimap2.smk"
include: "rules/DESeq2.smk"
include: "rules/volcano_plot.smk"

def get_contrast_names(wilds):
    subsets = comparisons.keys()
    outputs = []
    for subset in subsets:
        contrast_dir = checkpoints.dseq_dge.get(subset=subset).output[0]
        contrast_names = glob_wildcards(os.path.join(contrast_dir, "{sample}.txt"))
        outputs += expand(contrast_dir + "/{sample}.svg", sample=contrast_names.sample)
    return outputs

# Define rules that require minimal resources and can run on the same node
# where the actual snakemake is running. Should be low demand processess
localrules: run_all, dseq_initalize

# Rule run_all collects all outputs to force execution of the whole pipeline.
# Identified files will be produced
rule run_all:
    input:
        # alignments BAM/BAI files
        expand(
            SAMPLES_DIR + "/{s}/reads.{target}.sorted.{sufx}",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome'],
            sufx=['bam', 'bam.bai']),

        # flagstats
        expand(
            SAMPLES_DIR + "/{s}/reads.{target}.sorted.bam.flagstats.txt",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome']),

        # counts per transcript/gene
        expand(
            ANALYSIS_DIR +
            "/counts/{s}/reads.{target}.txt",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome'],
            ),
    
        # tables aggregated across tables
        expand(
            ANALYSIS_DIR + "/{subset}/all_genome_counts.txt",
            subset = comparisons.keys()
        ),
         # deseq2 differential gene expression analysis
        expand(
            ANALYSIS_DIR + "/{subset}/contrasts/permutations_list.txt",
            subset = comparisons.keys()
        ),        

        # expand(
        #     ANALYSIS_DIR + "/{subset}/contrasts/plots/{contrast}.svg",
        #     subset = comparisons.keys(),
        #     contrast = samples.keys()
        # ),

        get_contrast_names