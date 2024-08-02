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

# Set variables

SRA = config['SRA']

def count_lines(filename):
    with open(filename,'r') as file:
        line_count = sum(1 for line in file)//4
    return line_count

class SRA:
    def __init__(self, name, seq, group):
        self.name = name
        self.seq = seq
        self.group = group

Accessions = {}

for d in config['SRA']:
    s = SRA(*d)
    Accessions[s.name] = s

Ribo = []
RNA = []
for key,value in Accessions.items():
    if value.seq == "Ribo":
        Ribo.append(value.name)
    else:
        RNA.append(value.name)

# Define rules that require minimal resources and can run on the same node
# where the actual snakemake is running
localrules: run_all, aggregate_ATF4_reads, aggregate_cds_reads, count_lines_from_fastq, count_lines_from_trimmed_ribo

rule run_all:
    input: 
    # fetch fastq files
        expand(SAMPLES_DIR + "/{accession}.fastq",
        accession = Accessions.keys()),
    # wc file
        ANALYSIS_DIR + "/initial_count.txt",
    # trimmed riboseq reads
        expand(SAMPLES_DIR + "/{accession}/trimmed.fastq.gz",
        accession = Ribo),
    #rRNA free reads
        expand(EXP_DIR + "/{accession}/no_rrna.fastq",
        accession = Ribo),
    #rRNA reads
        expand(SAMPLES_DIR + "/{accession}/rrna.fastq.gz",
        accession = Ribo),
    #RiboTrim summary
        ANALYSIS_DIR + "/ribo_trim_count.txt",
    #Genome BAM
        expand(EXP_DIR + "/{accession}/Aligned.sortedByCoord.out.bam",
        accession = Accessions.keys()),
    #Genome Counts
        expand(EXP_DIR + "/{accession}/ReadsPerGene.out.tab",
        accession = Accessions.keys()),
    #Transcriptome BAM
        expand(EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.bam",
        accession = Accessions.keys()),
    #Transcriptome Sorted
        expand(EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam",
        accession = Accessions.keys()),
    #ATF4 PDF
        expand(EXP_DIR + "/{accession}/ATF4.pdf",
        accession = Accessions.keys()),
    #ATF4 Metaplot Txt
        expand(EXP_DIR + "/{accession}/ATF4.txt",
        accession = Accessions.keys()),  
    #Transcriptome BAI
        expand(EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam.bai",
        accession = Accessions.keys()),
    #ATF4 counts in a single tsv file
        ANALYSIS_DIR + "/ATF4_all_counts.txt",
    #Gene counts in a single tsv file
        ANALYSIS_DIR + "/all_sample_gene_counts.txt",
    #CDS counts in a single tsv file
        ANALYSIS_DIR + "/CDS_all_counts.txt",
    #CDS PDF
        expand(EXP_DIR + "/{accession}/CDS.pdf",
        accession = Accessions.keys()),   
    #CDS Metaplot Txt
        expand(EXP_DIR + "/{accession}/CDS.txt",
        accession = Accessions.keys()),  


rule fasterq_dump_from_SRA:
    output:
        SAMPLES_DIR + "/{accession}.fastq",
    log:
        LOGS_DIR + "/{accession}.log",
    params:
        accession = lambda wilds: wilds.accession,
    envmodules:
        "sratoolkit/3.0.10"    
    threads: 6
    resources:
        disk_mb = 50000,
        network = 1,
        mem_mb = 30000
    retries: 3
    shell:
        """
        fasterq-dump \
        {params.accession} \
        -e {threads} \
        -t /lscratch/$SLURM_JOBID \
        -o {output} \
        """

rule count_lines_from_fastq:
    input:
        fastqs = expand(SAMPLES_DIR + "/{accession}.fastq",accession = Accessions.keys()),
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
            pattern= r'(SRR\d+)\.fastq'
            match = re.search(pattern,file)
            dict[match.group(1)] = [match.group(1),
                    count_lines(file),
                    Accessions[match.group(1)].seq,
                    Accessions[match.group(1)].group
            ]
        df = pd.DataFrame(dict, index=None).T
        df.columns = ["Accession","lines","seq","group"]
        df.to_csv(output[0],sep='\t', index=False)

#Rule Cutadapt of RiboSeq footprints
rule Cutadapt_Ribo:
    input:
        SAMPLES_DIR + "/{accession}.fastq", 
    output:
        SAMPLES_DIR + "/{accession}/trimmed.fastq.gz",
    log:
        SAMPLES_DIR + "/{accession}/logs/cutadapt.log",
    envmodules:
        "cutadapt/4.4"
    shell:
        """
        cutadapt\
            -m 15 \
            -j 20 \
            -o {output} \
            {input} \
        """

#Rule Bowtie for filtering of rRNA reads
rule Bowtie2_for_rRNA_removal:
    input:
        fastq = SAMPLES_DIR + "/{accession}/trimmed.fastq.gz",
    output:
        mapped = SAMPLES_DIR + "/{accession}/rrna.fastq.gz",
        unmapped = EXP_DIR + "/{accession}/no_rrna.fastq",
    log:
        LOGS_DIR+ "{accession}.no_rrna.log"
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
            -L 20 \
            -x {params.genome} \
            --un {output.unmapped} \
            -U {input.fastq} \
            -S {output.mapped}
        """

rule count_lines_from_trimmed_ribo:
    input:
        unmapped = expand(EXP_DIR + "/{accession}/no_rrna.fastq",accession = Ribo),
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
            pattern = r'(SRR\d+)/'
            match = re.search(pattern,file)
            dict[match.group(1)] = [
                match.group(1),
                    count_lines(file),
                    Accessions[match.group(1)].seq,
                    Accessions[match.group(1)].group
            ]
        df = pd.DataFrame(dict, index=None).T
        df.columns = ["Accession","lines","seq","group"]
        df.to_csv(output[0],sep='\t', index=False)

def input_for_alignment_rules(accession):
    if accession.seq == "Ribo": 
        return os.path.join(
            EXP_DIR, accession.name, "no_rrna.fastq")

    return os.path.join(
        SAMPLES_DIR, accession.name + ".fastq")


#Rule Align RiboSeq with STAR to Genome and Transcriptome
rule STAR_Riboseq:
    input:
        fastq = lambda ws: input_for_alignment_rules(Accessions[ws.accession]),
        genome = DATA_DIR + "/STAR_Index_108"
    output:
        Genome_BAM = EXP_DIR + "/{accession}/Aligned.sortedByCoord.out.bam",
        Transcriptome_SAM = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.bam",
        COUNTS = EXP_DIR + "/{accession}/ReadsPerGene.out.tab",
    log:
        LOGS_DIR + "/{accession}.STAR.log"
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
            --outFileNamePrefix experiment/{wildcards.accession}/ \
            --quantMode GeneCounts TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 31532137230 \
            --outSAMattributes All
        """

## For ribo-seq, reads are unstranded (column 1 of ReadsPerGene)
rule aggregate_gene_counts_ribo:
    input:
        counts = expand(EXP_DIR + "/{accession}/ReadsPerGene.out.tab", accession = Ribo)
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
        counts = expand(EXP_DIR + "/{accession}/ReadsPerGene.out.tab", accession=RNA)
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
        EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.bam",
    output:
        EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam",
    envmodules:
        "samtools/1.19"
    resources: 
        runtime=100, 
        mem_mb=20000, 
        disk_mb=20000,
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
        EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam",
    output:
        EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam.bai",
    envmodules:
        "samtools/1.19"
    resources: 
        runtime=100, 
        mem_mb=20000, 
        disk_mb=20000,
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

# count_aligned_reads_per_transcript counts the reads aligned on each transcript.
rule count_aligned_reads_per_transcript:
    input:
        aligned = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.bam",
        transcript_tab = DATA_DIR + "/hg38/transcript-gene.tab",
    output:
        EXP_DIR + "/{accession}/ReadsPerTranscript.tab",
    conda:
        "envs/bam.yml" 
    shell:
        """

        {SCRIPTS_DIR}/sam_per_ref_count_statistics.py \
            --ifile {input.aligned} \
            --ref-col-name transcript \
            --cnt-col-name count \
            --opt-col-name sample \
            --opt-col-val {wildcards.accession} \
            | table-join.py \
                --table1 - \
                --table2 {input.transcript_tab} \
                --key1 transcript \
                --key2 transcript \
                > {output}
        """

#Metaplot of ATF4
rule ATF4_Metaplot:
    input:
        BAM = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam",
        BAI = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam.bai",
    output:
        PDF = EXP_DIR + "/{accession}/ATF4.pdf",
        TXT = EXP_DIR + "/{accession}/ATF4.txt",
    params:
        BED = DATA_DIR + "/ATF4_mRNA.bed"
    log:
        LOGS_DIR + "/{accession}.metaplot.log"
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
         expand(EXP_DIR + "/{accession}/ATF4.txt", accession = Ribo),
    output:
        ANALYSIS_DIR + "/ATF4_all_counts.txt",
    run:
        import pandas as pd
        import re
        dfs= []
        for file_path in input:
           pattern = r'(SRR\d+)/'
           name = re.search(pattern,file_path)
           df = pd.read_csv(file_path, delimiter='\t')
           df.rename(columns={'count': name.group(1)}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)

#Metaplot of All genes
rule CDS_Metaplot:
    input:
        BAM = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam",
        BAI = EXP_DIR + "/{accession}/Aligned.toTranscriptome.out.sorted.bam.bai",
    output:
        PDF = EXP_DIR + "/{accession}/CDS.pdf",
        TXT = EXP_DIR + "/{accession}/CDS.txt",
    params:
        BED = DATA_DIR + "/cds_genic.bed"
    log:
        LOGS_DIR + "/{accession}.cds.metaplot.log"
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
         expand(EXP_DIR + "/{accession}/CDS.txt", accession = Ribo),
    output:
        ANALYSIS_DIR + "/CDS_all_counts.txt",
    run:
        import pandas as pd
        import re
        dfs= []
        for file_path in input:
           pattern = r'(SRR\d+)/'
           name = re.search(pattern,file_path)
           df = pd.read_csv(file_path, delimiter='\t')
           df.rename(columns={'count': name.group(1)}, inplace=True)
           dfs.append(df)
        combined_df = dfs[0]
        for df in dfs[1:]:
            combined_df = pd.merge(combined_df , df, on = "pos", how = "inner")
        combined_df.to_csv(output[0], sep = '\t', index=False)
