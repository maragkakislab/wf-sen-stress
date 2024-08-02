import sys
import glob

config_path = "config/config.yml"
configfile: config_path



# Set directories
SAMPLES_DIR = config['SAMPLES_DIR']
ANALYSIS_DIR = config['ANALYSIS_DIR']
LOGS_DIR = config['LOGS_DIR']
DATA_DIR = config['DATA_DIR']
EXP_DIR = config['EXP_DIR']
SCRIPTS_DIR = config['SCRIPTS_DIR']
RAW_DIR = config['RAW_DIR']

# Set variables


def count_lines(filename):
    with open(filename,'r') as file:
        line_count = sum(1 for line in file)//4
    return line_count

class RSeq:
    def __init__(self, name, seq, group):
        self.name = name
        self.seq = seq
        self.group = group

samples = {}

for d in config['RAW']:
    s = RSeq(*d)
    samples[s.name] = s

Ribo = []
RNA = []
for key,value in samples.items():
    if value.seq == "Ribo":
        Ribo.append(value.name)
    else:
        RNA.append(value.name)


rule run_all:
    input: 
    # wc file
        ANALYSIS_DIR + "/initial_count.txt",
    # trimmed riboseq reads
        expand(SAMPLES_DIR + "/{sample}/trimmed.fastq.gz",
        sample = Ribo),
    #rRNA free reads
        expand(EXP_DIR + "/{sample}/no_rrna.fastq",
        sample = Ribo),
    #rRNA reads
        expand(SAMPLES_DIR + "/{sample}/rrna.fastq.gz",
        sample = Ribo),
    #RiboTrim summary
        ANALYSIS_DIR + "/ribo_trim_count.txt",
    #Genome BAM
        expand(EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam",
        sample = samples.keys()),
    #Genome Counts
        expand(EXP_DIR + "/{sample}/ReadsPerGene.out.tab",
        sample = samples.keys()),
    #Transcriptome BAM
        expand(EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.bam",
        sample = samples.keys()),
    #Transcriptome Sorted
        expand(EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam",
        sample = samples.keys()),
    # #ATF4 Clean PDF
    #     expand(EXP_DIR + "/{sample}/ATF4_min.pdf",
    #     sample = samples.keys()),
    # #ATF4 Clean Metaplot Txt
    #     expand(EXP_DIR + "/{sample}/ATF4_min.txt",
    #     sample = samples.keys()),
    #ATF4 PDF
        expand(EXP_DIR + "/{sample}/ATF4.pdf",
        sample = samples.keys()),
    #Feature count summary
        ANALYSIS_DIR + "/count_features.txt",
    #ATF4 Metaplot Txt
        expand(EXP_DIR + "/{sample}/ATF4.txt",
        sample = samples.keys()),
    #Riboseq Quality Control
        expand(EXP_DIR + "/{sample}/Aligned.sortedByCoord.out_qual.pdf",
        sample = Ribo),               
         
    #ATF4 Ribo counts in a single tsv file
        ANALYSIS_DIR + "/ATF4_all_counts.txt",
    # #ATF4 Ribo counts in a single tsv file   
    #     ANALYSIS_DIR + "/ATF4_all_clean_counts.txt",
    #Gene counts in a single tsv file
        ANALYSIS_DIR + "/all_sample_gene_counts.txt",
    #CDS counts in a single tsv file
        ANALYSIS_DIR + "/CDS_all_counts.txt",
    # #CDS PDF
    #     expand(EXP_DIR + "/{sample}/CDS_clean.pdf",
    #     sample = samples.keys()),   
    # #CDS Metaplot Txt
    #     expand(EXP_DIR + "/{sample}/CDS_clean.txt",
    #     sample = samples.keys()),  
    # #No Dups ATF4
    #     expand(EXP_DIR + "/{sample}/ATF4_nodups.pdf",
    #     sample = samples.keys()),   
    # #No Dups ATF4
    #     expand(EXP_DIR + "/{sample}/ATF4_nodups.txt",
    #     sample = samples.keys()),  
    # #ATF4 No Dups Ribo counts in a single tsv file   
    #     ANALYSIS_DIR + "/ATF4_all_nodups_counts.txt",




rule count_lines_from_fastq:
    input:
        fastqs = expand(RAW_DIR + "/{sample}.fastq",sample = samples.keys()),
    output:
        ANALYSIS_DIR + "/initial_count.txt",
    resources:
        disk_mb = 5000,
        mem_mb = 5000,
    run:
        import re
        import pandas as pd
        files = input.fastqs   
        dict = {}
        for file in files:
            pattern= r'/([^/]+)\.fastq'
            match = re.search(pattern,file)
            dict[match.group(1)] = [match.group(1),
                    count_lines(file),
                    samples[match.group(1)].seq,
                    samples[match.group(1)].group
            ]
        df = pd.DataFrame(dict, index=None).T
        df.columns = ["sample","lines","seq","group"]
        df.to_csv(output[0],sep='\t')

#Rule Cutadapt of RiboSeq footprints
rule Cutadapt_Ribo:
    input:
        RAW_DIR + "/{sample}.fastq", 
    output:
        SAMPLES_DIR + "/{sample}/trimmed.fastq.gz",
    log:
        SAMPLES_DIR + "/{sample}/logs/cutadapt.log",
    envmodules:
        "cutadapt/4.4"
    params:
        adaptor = config['ADAPTOR']
    threads:
        8
    shell:
        """
        cutadapt \
            -m 22 \
            -M 34 \
            -a {params.adaptor} \
            -j {threads} \
            -o {output} \
            {input} \
        """

#Rule Bowtie for filtering of rRNA reads
rule Bowtie2_for_rRNA_removal:
    input:
        fastq = SAMPLES_DIR + "/{sample}/trimmed.fastq.gz",
    output:
        mapped = SAMPLES_DIR + "/{sample}/rrna.fastq.gz",
        unmapped = EXP_DIR + "/{sample}/no_rrna.fastq",
    log:
        LOGS_DIR+ "{sample}.no_rrna.log"
    threads: 
        12
    envmodules:
        "bowtie/2-2.5.1"
    params:
        genome = DATA_DIR + "/Index_Bowtie/rRNA/rRNA",
    threads: 50
    resources:
        mem_mb = 100*1024,
        runtime = 2*24*60
    shell:
        """
        bowtie2\
            -p {threads}\
            -x {params.genome} \
            --un {output.unmapped} \
            -U {input.fastq} \
            -S {output.mapped}
        """

rule count_lines_from_trimmed_ribo:
    input:
        unmapped = expand(EXP_DIR + "/{sample}/no_rrna.fastq",sample = Ribo),
    output:
        ANALYSIS_DIR + "/ribo_trim_count.txt",
    resources:
        disk_mb = 5000,
        mem_mb = 5000,
    run:
        import re
        import pandas as pd
        files = input.unmapped  
        dict = {}
        for file in files:
            pattern= r'/([^/]+)/no_rrna\.fastq'
            match = re.search(pattern,file)
            dict[match.group(1)] = [match.group(1),
                    count_lines(file),
                    samples[match.group(1)].seq,
                    samples[match.group(1)].group
            ]
        df = pd.DataFrame(dict, index=None).T
        df.columns = ["sample","lines","seq","group"]
        df.to_csv(output[0],sep='\t')

def input_for_alignment_rules(sample):
    if sample.seq == "Ribo": 
        return os.path.join(
            EXP_DIR, sample.name, "no_rrna.fastq")

    return os.path.join(
        RAW_DIR, sample.name + ".fastq")

#Rule Align RiboSeq with STAR to Genome and Transcriptome
rule STAR_Riboseq:
    input:
        fastq = lambda ws: input_for_alignment_rules(samples[ws.sample]),
        genome = DATA_DIR + "/STAR_Index_108"
    output:
        Genome_BAM = EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam",
        Transcriptome_SAM = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.bam",
        COUNTS = EXP_DIR + "/{sample}/ReadsPerGene.out.tab",
    log:
        LOGS_DIR + "/{sample}.STAR.log"
    threads: 
        12
    envmodules:
        "STAR/2.7.10b"
    resources: 
        runtime=100, 
        mem_mb=100*1024, 
        disk_mb=100*1024,
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 51 \
            --outFilterMismatchNmax 2 \
            --genomeDir {input.genome} \
            --readFilesIn {input.fastq} \
            --outFileNamePrefix experiment/{wildcards.sample}/ \
            --quantMode GeneCounts TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 31532137230 \
            --outSAMattributes All
        """

#FeatureCounts
rule FeatureCounts_CDS:
    input: 
        bam = EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam",
        gtf = DATA_DIR + "/hg38/Homo_sapiens.GRCh38.111.gtf"
    output:
        counts = EXP_DIR + "/{sample}/counts.cds.txt",
        summary = EXP_DIR + "/{sample}/counts.cds.txt.summary",
    envmodules:
        "subread/2.0.3"
    threads:
        5
    shell:
        """
        featureCounts \
        -T {threads} \
        -t CDS \
        -d 15 \
        -g gene_id \
        -a {input.gtf} \
        -o {output.counts} \
        -M \
        {input.bam}
        """

#FeatureCounts
rule FeatureCounts_utr:
    input: 
        bam = EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam",
        gtf = DATA_DIR + "/hg38/Homo_sapiens.GRCh38.111.gtf"
    output:
        counts = EXP_DIR + "/{sample}/counts.utr.txt",
        summary = EXP_DIR + "/{sample}/counts.utr.txt.summary",
    envmodules:
        "subread/2.0.3"
    threads:
        5
    shell:
        """
        featureCounts \
        -T {threads} \
        -t three_prime_utr,five_prime_utr \
        -d 15 \
        -g gene_id \
        -a {input.gtf} \
        -o {output.counts} \
        {input.bam}
        """

#Create count summary table
rule count_summary:
    input:
        cds = expand(EXP_DIR + "/{sample}/counts.cds.txt.summary", sample = samples.keys()),
        utr = expand(EXP_DIR + "/{sample}/counts.utr.txt.summary", sample = samples.keys()),
    output:
        ANALYSIS_DIR + "/count_features.txt",
    run:
        import re
        import pandas as pd
        cds = input.cds
        utr = input.utr
        dict_cds = {}
        dict_utr = {}
        for file_path in cds:
            with open(file_path, 'r') as file:
                name = file_path.split('/')[-2]
                lines = file.readlines()
                parts = lines[1].split('\t')
                aligned = int(parts[1])
                dict_cds[name] = aligned
        for file_path in utr:
            with open(file_path, 'r') as file:
                name = file_path.split('/')[-2]
                lines = file.readlines()
                parts = lines[1].split('\t')
                aligned = int(parts[1])
                dict_utr[name] = aligned
        series_utr = pd.Series(dict_utr, name='dict_utr')
        series_cds = pd.Series(dict_cds, name='dict_cds') 
        df = pd.concat([series_utr, series_cds], axis=1)           
        df.to_csv(output[0], sep = '\t')  
                





#RiboSeq Quality Control with Ribotish
rule ribotish:
    input: 
        bam = EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam",
        bai = EXP_DIR + "/{sample}/Aligned.sortedByCoord.out.bam.bai",
        gtf = DATA_DIR + "/hg38/genes.gtf",
    output: 
        EXP_DIR + "/{sample}/Aligned.sortedByCoord.out_qual.pdf"
    conda:
        "envs/ribotish.yml" 
    resources: 
        runtime=100, 
        mem_mb=100*1024, 
        disk_mb=100*1024,    
    shell:
        """
        ribotish quality \
            -b {input.bam} \
            -g {input.gtf} \
        """

## For ribo-seq, reads are unstranded (column 1 of ReadsPerGene)
rule aggregate_gene_counts_ribo:
    input:
        counts = expand(EXP_DIR + "/{sample}/ReadsPerGene.out.tab", sample = Ribo)
    output:
        ANALYSIS_DIR + "/RNA_sample_gene_counts.txt"
    run:
        import pandas as pd
        import re
        files = input.counts  
        dict = {}
        dfs= []
        for file in files:
            pattern= r'/([^/]+)/ReadsPerGene\.out\.tab'
            match = re.search(pattern,file)       
            df = pd.read_csv(file, delimiter='\t')
            df = df.drop(range(4)).drop(df.columns[[2,3]], axis=1).reset_index(drop=True)
            df.columns = ["ENSG", match.group(1)]
            dfs.append(df)
        combined_df = dfs[0]      
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df, df, on = "ENSG", how = "inner")   
        combined_df.reset_index(drop=True, inplace=True) 
        combined_df.to_csv(output[0], sep = '\t', index=False)            

## For rna-seq, reads are stranded reverse (column 3 of ReadsPerGene)       
rule aggregate_gene_counts_RNA:
    input:
        counts = expand(EXP_DIR + "/{sample}/ReadsPerGene.out.tab", sample=RNA)
    output:
        ANALYSIS_DIR + "/ribo_sample_gene_counts.txt"
    run:
        import pandas as pd
        import re
        files = input.counts  
        dict = {}
        dfs= []
        for file in files:
            pattern= r'/([^/]+)/ReadsPerGene\.out\.tab'
            match = re.search(pattern,file)       
            df = pd.read_csv(file, delimiter='\t')
            df = df.drop(range(4)).drop(df.columns[[1,2]], axis=1).reset_index(drop=True)
            df.columns = ["ENSG", match.group(1)]
            dfs.append(df)
        combined_df = dfs[0]      
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df, df, on = "ENSG", how = "inner")   
        combined_df.reset_index(drop=True, inplace=True) 
        combined_df.to_csv(output[0], sep = '\t', index=False)            


rule aggregate_all_counts:
    input: 
        Ribo_counts = ANALYSIS_DIR + "/ribo_sample_gene_counts.txt",
        RNA_counts = ANALYSIS_DIR + "/RNA_sample_gene_counts.txt",
    output:
        ANALYSIS_DIR + "/all_sample_gene_counts.txt",
    run:
        import pandas as pd
        Ribo_counts = pd.read_csv(input.Ribo_counts, delimiter='\t')
        RNA_counts = pd.read_csv(input.RNA_counts, delimiter='\t')
        all_counts = pd.merge(Ribo_counts, RNA_counts, on="ENSG",how='inner')
        all_counts.to_csv(output[0], sep = '\t', index=False)  
        
#Sort Bam files
rule bam_sort:
    input:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.bam",
    output:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam",
    envmodules:
        "samtools/1.19"
    resources: 
        runtime=100, 
        mem_mb=50000, 
        disk_mb=50000,
    threads:
        4
    shell:
        """
        samtools \
            sort \
            -@ {threads} \
            -o {output} \
            {input}

        """

#Index Bam files
rule bam_index:
    input:
        EXP_DIR + "/{sample}/{prefix}.bam",
    output:
        EXP_DIR + "/{sample}/{prefix}.bam.bai",
    envmodules:
        "samtools/1.19"
    resources: 
        runtime=100, 
        mem_mb=50000, 
        disk_mb=50000,
    threads:
        4
    shell:
        """
        samtools \
            index \
            -b \
            -@ {threads} \
            -o {output} \
            {input}
        """

# count_aligned_reads_per_transcript counts the reads aligned on each
# transcript.
rule count_aligned_reads_per_transcript:
    input:
        aligned = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.bam",
        transcript_tab = DATA_DIR + "/hg38/transcript-gene.tab",
    output:
        EXP_DIR + "/{sample}/ReadsPerTranscript.tab",
    conda:
        "envs/bam.yml" 
    shell:
        """

        {SCRIPTS_DIR}/sam_per_ref_count_statistics.py \
            --ifile {input.aligned} \
            --ref-col-name transcript \
            --cnt-col-name count \
            --opt-col-name sample \
            --opt-col-val {wildcards.sample} \
            | table-join.py \
                --table1 - \
                --table2 {input.transcript_tab} \
                --key1 transcript \
                --key2 transcript \
                > {output}
        """


#Mark putative PCR duplicates using Picard
rule duplicates_cleanup:
    input:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam",
    output:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.nodups.bam",
    envmodules:
        "samtools/1.19"
    resources: 
        runtime=100, 
        mem_mb=50000, 
        disk_mb=50000,
    threads:
        4
    shell:
        """
        samtools \
            view \
            -F 256 \
            -@ {threads} \
            -o {output} \
            {input}
        """

rule bam_min_length:
    input:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam",
    output:
        EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.minlen.bam",
    resources: 
        runtime=100, 
        mem_mb=50000, 
        disk_mb=50000,
    conda:
        "envs/bam.yml" 
    shell:
        """

        {SCRIPTS_DIR}/minlen.py \
            --ifile {input} \
            --out {output} \
            --min 28

        """


#Metaplot of ATF4
rule ATF4_clean_Metaplot:
    input:
        BAM = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.minlen.bam",
        BAI = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.minlen.bam.bai",
    output:
        PDF = EXP_DIR + "/{sample}/ATF4_min.pdf",
        TXT = EXP_DIR + "/{sample}/ATF4_min.txt",
    params:
        BED = DATA_DIR + "/ATF4_mRNA.bed"
    log:
        LOGS_DIR + "/{sample}.metaplot.log"
    conda:
        "envs/biometatools.yml"
    shell:
        """
         rel-pos-bam-vs-bed.py \
            --bam {input.BAM} \
            --bed {params.BED} \
            --pdf {output.PDF} \
            --up 100 \
            --down 1500 \
            > {output.TXT}
        """


####Aggregate ATF4 metaplot data
rule aggregate_clean_ATF4_reads:
    input:
         expand(EXP_DIR + "/{sample}/ATF4_min.txt", sample = Ribo),
    output:
        ANALYSIS_DIR + "/ATF4_all_clean_counts.txt",
    run:
        import pandas as pd
        import re
        dfs= []
        for file_path in input:
           name = file_path.split('/')[-2]
           df = pd.read_csv(file_path, delimiter='\t')
           df.rename(columns={'count': name}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)



#Metaplot of ATF4
rule ATF4_Metaplot:
    input:
        BAM = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam",
        BAI = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.bam.bai",
    output:
        PDF = EXP_DIR + "/{sample}/ATF4.pdf",
        TXT = EXP_DIR + "/{sample}/ATF4.txt",
    params:
        BED = DATA_DIR + "/ATF4_mRNA.bed"
    log:
        LOGS_DIR + "/{sample}.metaplot.log"
    conda:
        "envs/biometatools.yml"
    shell:
        """
         rel-pos-bam-vs-bed.py \
            --bam {input.BAM} \
            --bed {params.BED} \
            --pdf {output.PDF} \
            --up 100 \
            --down 1500 \
            > {output.TXT}
        """


####Aggregate ATF4 metaplot data
rule aggregate_ATF4_reads:
    input:
         expand(EXP_DIR + "/{sample}/ATF4.txt", sample = Ribo),
    output:
        ANALYSIS_DIR + "/ATF4_all_counts.txt",
    run:
        import pandas as pd
        import re
        dfs= []
        for file_path in input:
           name = file_path.split('/')[-2]
           df = pd.read_csv(file_path, delimiter='\t')
           df.rename(columns={'count': name}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)


#Metaplot of All genes
rule CDS_Metaplot:
    input:
        BAM = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.nodups.bam",
        BAI = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.nodups.bam.bai",
    output:
        PDF = EXP_DIR + "/{sample}/CDS_clean.pdf",
        TXT = EXP_DIR + "/{sample}/CDS_clean.txt",
    params:
        BED = DATA_DIR + "/cds_genic.bed"
    log:
        LOGS_DIR + "/{sample}.cds.metaplot.log"
    conda:
        "envs/biometatools.yml"
    resources: 
        runtime=100, 
        mem_mb=40000, 
        disk_mb=40000,
    shell:
        """
         rel-pos-bam-vs-bed.py \
            --bam {input.BAM} \
            --bed {params.BED} \
            --pdf {output.PDF} \
            --up 50 \
            --down 500 \
            > {output.TXT}
        """

####Aggregate cds metaplot data
rule aggregate_cds_reads:
    input:
         fastq = expand(EXP_DIR + "/{sample}/CDS.txt", sample = Ribo),
    output:
        ANALYSIS_DIR + "/CDS_all_counts.txt",
    run:
        import pandas as pd
        import re
        files = input.fastq
        dfs= []
        for file in files:
           name = file.split('/')[-2]
           df = pd.read_csv(file, delimiter='\t')
           df.rename(columns={'count': name}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)


#Metaplot of ATF4
rule ATF4_nodups_Metaplot:
    input:
        BAM = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.nodups.bam",
        BAI = EXP_DIR + "/{sample}/Aligned.toTranscriptome.out.sorted.nodups.bam.bai",
    output:
        PDF = EXP_DIR + "/{sample}/ATF4_nodups.pdf",
        TXT = EXP_DIR + "/{sample}/ATF4_nodups.txt",
    params:
        BED = DATA_DIR + "/ATF4_mRNA.bed"
    log:
        LOGS_DIR + "/{sample}.metaplot.log"
    conda:
        "envs/biometatools.yml"
    shell:
        """
         rel-pos-bam-vs-bed.py \
            --bam {input.BAM} \
            --bed {params.BED} \
            --pdf {output.PDF} \
            --up 100 \
            --down 1500 \
            > {output.TXT}
        """


####Aggregate ATF4 metaplot data
rule aggregate_nodups_ATF4_reads:
    input:
         expand(EXP_DIR + "/{sample}/ATF4_nodups.txt", sample = Ribo),
    output:
        ANALYSIS_DIR + "/ATF4_all_nodups_counts.txt",
    run:
        import pandas as pd
        import re
        dfs= []
        for file_path in input:
           name = file_path.split('/')[-2]
           df = pd.read_csv(file_path, delimiter='\t')
           df.rename(columns={'count': name}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)
