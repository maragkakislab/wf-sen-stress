
def deseq_counts(subset):
    samples = comparisons[subset]['samples']
    files = [ANALYSIS_DIR + f"/counts/{sample[0]}/reads.toGenome.txt" for sample in samples]
    return files

rule dseq_initalize:
    input:
        files = lambda wilds: deseq_counts(wilds.subset),
    output:
        counts = ANALYSIS_DIR + "/{subset}/all_genome_counts.txt",
        metafile = ANALYSIS_DIR + "/{subset}/metafile.txt",
    run:
        # create a metafile
        subset = wildcards.subset
        print(f"Subset: {subset}")
        samples = comparisons[subset]['samples']
        cols = comparisons[subset]['cols']
        print(f"cols: {cols}")
        data = pd.DataFrame(columns = cols)
        for i in range(len(samples)):
            data.loc[i] = samples[i]
        data.to_csv(output.metafile, sep = '\t',index=False)
       
        # aggregate reads
        dfs = [pd.read_csv(file_path, delimiter='\t') for file_path in input['files']]
        combined_df = pd.concat(dfs, axis=0, ignore_index=True)
        wide_df = combined_df.pivot(index= "gene", columns= 'sample', values = 'count')
        wide_df.to_csv(output.counts, sep = '\t')


rule ensembl_to_geneID:
    input:
        transcript_tab =  DATA_DIR + "/" + ASSEMBLY + "/transcript_gene.tab",
    output:
        ENSG_metadata =  DATA_DIR + "/" + ASSEMBLY + "/ENSG_metadata.tab",
    resources:
        mem_mb=5*1024,
        runtime = 1*60
    retries:
        3
    params:
        identifier = "hgnc_symbol",  # this will be moved to the config subset dictionary
        genome = "hsapiens_gene_ensembl", # this will be moved to the config subset dictionary
    envmodules:
        "R/4.4.0",
    script:
        "../scripts/biomart.R"

#localrules: dseq_dge ###toggle to debug this in dry run
checkpoint dseq_dge:
    input:
        counts = ANALYSIS_DIR + "/{subset}/all_genome_counts.txt",
        metafile = ANALYSIS_DIR + "/{subset}/metafile.txt",
        ENSG_metadata =  DATA_DIR + "/" + ASSEMBLY + "/ENSG_metadata.tab",
    output:
       out = directory(ANALYSIS_DIR + "/{subset}/contrasts/"),
       permutations = ANALYSIS_DIR + "/{subset}/permutations_list.txt",
    resources:
        mem_mb=5*1024,
        runtime=3*60,
    params:
        model = " ~ condition",
        delim = "\t",
        odir = ANALYSIS_DIR + "/{subset}/contrasts/",
    envmodules:
        "R/4.4.0",
    script:
        "../scripts/deseq2.R"
