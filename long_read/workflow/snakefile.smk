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
SAMPLES_DIR = config['SAMPLES_DIR']
SCRIPTS_DIR = config['SCRIPTS_DIR']
RSYNC_PATH = config['RSYNC_PATH']



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


##### include rules #####
include: "rules/dorado.smk"
include: "rules/minimap2.smk"

# Define rules that require minimal resources and can run on the same node
# where the actual snakemake is running. Should be low demand processess
localrules: run_all, aggregate

# Rule run_all collects all outputs to force execution of the whole pipeline.
# Identified files will be produced
rule run_all:
    input:
        # Alignments BAM/BAI files
        expand(
            SAMPLES_DIR + "/{s}/reads.{target}.sorted.{sufx}",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome'],
            sufx=['bam', 'bam.bai']),

        # Flagstats
        expand(
            SAMPLES_DIR + "/{s}/reads.{target}.sorted.bam.flagstats.txt",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome']),

        # Counts per transcript/gene
        expand(
            ANALYSIS_DIR +
            "/counts/{s}/reads.{target}.txt",
            s=samples.keys(),
            target=['toGenome', 'toTranscriptome'],
            ),
 
        # tables aggregated across tables
        expand(ANALYSIS_DIR + "/{e}/all_genome_counts.txt",
        e = experiments.keys())
        ####ANALYSIS_DIR + "/all_transcriptome_counts.txt",

############## Rules

# ###### Generate a single dataframe

rule aggregate:
    input:
       file_paths = expand(ANALYSIS_DIR + "/counts/{s}/reads.toGenome.txt", s=samples.keys()),
    output:
        expand(ANALYSIS_DIR + "/{e}/all_genome_counts.txt",
        e=experiments.keys()),
    run:
        import pandas as pd
        #### makes a dictionary where experiments are keys 
        #### and file_paths to counts are items 
        expts = {key:[] for key in experiments.keys()}
        for file_path in input['file_paths']:
            print(file_path)
            split = file_path.split('/')
            if len(split) > 1:
                s = split[-2]
                print(s)
                e = samples[s].parent_exp 
                if e in expts.keys():
                    expts[e].append(file_path)

        #### loops aggretation over each experiment and saves it to an experiment folder
        for expt in expts.keys():
            dfs = [pd.read_csv(file_path, delimiter='\t') for file_path in expts[expt]]
            combined_df = pd.concat(dfs, axis=0, ignore_index=True)
            wide_df = combined_df.pivot(index= "gene", columns= 'sample', values = 'count')
            wide_df.to_csv(ANALYSIS_DIR + f"/{expt}_all_genome_counts.txt", sep = '\t')


# rule transcriptome_combined:
#     input:
#          expand(ANALYSIS_DIR + "/counts/{s}/reads.toTranscriptome.txt", s=samples.keys()),
#     output:
#         ANALYSIS_DIR + "/all_transcriptome_counts.txt",
#     run:
#         import pandas as pd
#         dfs = [pd.read_csv(file_path, delimiter='\t') for file_path in input]
#         combined_df = pd.concat(dfs, axis=0, ignore_index=True)
#         wide_df = combined_df.pivot(index= "transcript", columns= 'sample', values = 'count')
#         wide_df.to_csv(output[0], sep = '\t')





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

