import os
import glob

SAMPLES = sorted(list(set([os.path.basename(f).split('_')[0] for f in glob.glob("raw_data/*_1.fastq.gz")])))

rule all:
    input:
         expand("results/fastp/{sample}_fastp.json", sample=SAMPLES),
         expand("results/bracken/{sample]_bracken.tsv", sample=SAMPLES),
         expand("results/bracken/{sample]_cazy.tsv", sample=SAMPLES),
         "results/taxa_coda_final.csv",
         "results/cazy_coda_final.csv",
         "results/target_gh_final.csv",
         "results/taxa_coda_final.csv",
        "results/cazy_coda_final.csv",
        "results/target_gh_final.csv",
        "results/taxon_cazy_integration.csv",
        "results/cazy_aitchison_pca.png",
        "results/Fig1_Taxa_Volcano.png",
        "results/Fig2_Top20_Taxa_Barplot.png",
        "results/Fig4_Target_GH_Barplot.png",
        "results/Fig5_Taxon_CAZy_Heatmap.png"

rule fastp:
    input:
         r1 = "raw_data/{sample}_1.fastq.gz",
         r2 = "raw_data/{sample}_2.fastq.gz"
    output:
         report = "results/kraken2/{sample}_kraken2.report",
         output = temp("results/kraken2/{sample}_kraken2.output")
         threads: 8
    params:
        db = "databases/kraken2_std_8gb"
    shell:
        "kraken2 --db {params.db} --threads {threads} --report {output.report} --paired {input.r1} {input.r2} > {output.output}"

rule bracken:
      input:
        report = "results/kraken2/{sample}_kraken2.report"
    output:
        abundances = "results/bracken/{sample}_bracken.tsv"
    threads: 1
    params:
        db = "databases/kraken2_std_8gb",
        level = "S"
    shell:
        "bracken -d {params.db} -i {input.report} -o {output.abundances} -r 150 -l {params.level}"

rule diamond:
    input:
        r1 = "results/qc/{sample}_1.trimmed.fastq.gz",
        r2 = "results/qc/{sample}_2.trimmed.fastq.gz"
    output:
        tsv = "results/diamond/{sample}_cazy.tsv"
    threads: 8
    params:
        db = "databases/cazy_diamond/CAZyDB"
    shell:
        """
        zcat {input.r1} {input.r2} | diamond blastx -d {params.db} -q /dev/stdin \
        -o {output.tsv} -p {threads} --block-size 1 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue \
        --max-target-seqs 1 --evalue 1e-10 --id 60 --query-cover 80
        """

rule coda_statistical_modeling:
    input:
        tax = "results/merged_bracken_species.csv",
        cazy = "results/merged_cazy_rpk.csv",
        meta = "metadata.csv"
    output:
        "results/taxa_coda_final.csv",
        "results/cazy_coda_final.csv",
        "results/target_gh_final.csv",
        "results/taxon_cazy_integration.csv",
        "results/cazy_aitchison_pca.png"
    shell:
        "Rscript master_coda_pipeline.R"

rule generate_publication_figures:
    input:
        "results/taxa_coda_final.csv",
        "results/target_gh_final.csv",
        "results/taxon_cazy_integration.csv"
    output:
        "results/Fig1_Taxa_Volcano.png",
        "results/Fig2_Top20_Taxa_Barplot.png",
        "results/Fig4_Target_GH_Barplot.png",
        "results/Fig5_Taxon_CAZy_Heatmap.png"
    shell:
        "python3 generate_final_figures.py"

         
