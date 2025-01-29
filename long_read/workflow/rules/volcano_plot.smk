
    
# def get_contrast_names(wilds):
#     subsets = comparisons.keys()
#     contrast_dir = expand(checkpoints.get_contrast_names.get(**wilds).output[0], subset = subsets)
#     contrast_names = glob_wildcards(os.path.join(contrast_dir, "{sample}.txt"))
#     return expand(ANALYSIS_DIR + "/{subset}/contrasts/plots/{SAMPLE}.svg", 
#     SAMPLE=contrast_names, subset = subsets)

rule volcano_plot:
    input:
        tsv_data = ANALYSIS_DIR + "/{subset}/contrasts/{contrast}.txt",
    output:
        plots = ANALYSIS_DIR + "/{subset}/contrasts/{contrast}.svg",
    resources:
        mem_mb=5*1024,
    params:
        odir = ANALYSIS_DIR + "/{subset}/contrasts/{contrast}/",
        baseline = 1,
        top_number = 10, #Move to config
    envmodules:
        "R/4.4.0",
    script:
        "../scripts/volcano_plot.R"