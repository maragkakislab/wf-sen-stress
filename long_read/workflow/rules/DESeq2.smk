
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


rule dseq_dge:
    input:
        counts = ANALYSIS_DIR + "/{subset}/all_genome_counts.txt",
        metafile = ANALYSIS_DIR + "/{subset}/metafile.txt",
    output:
       directory(ANALYSIS_DIR + "/{subset}/contrasts/")
    params:
        model = " ~ condition",
        delim = "\t",
        odir = ANALYSIS_DIR + "/{subset}/contrasts/"
    envmodules:
        "R/4.4.0",
    script:
        "../scripts/deseq2.R"