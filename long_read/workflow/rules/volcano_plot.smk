
    

rule volcano_plot:
    input:
        tsv_data = ANALYSIS_DIR + "/{subset}/contrasts/{contrast}.dge.txt",
    output:
        plot = ANALYSIS_DIR + "/{subset}/contrasts/plots/{contrast}/volcano.svg",
    resources:
        mem_mb=5*1024,
    wildcard_constraints:  
        contrast = "[^/]+", # This prevents {contrast} from becoming a directory path
    params:
        odir = ANALYSIS_DIR + "/{subset}/contrasts/plots/{contrast}/",
        baseline = 1,
        top_number = 10, #Move to config
        contrast = lambda wilds: wilds.contrast
    envmodules:
        "R/4.4.0",
    script:
        "../scripts/volcano_plot.R"