from snakemake.utils import min_version
min_version("6.0")
configfile: "resources/config.yaml"
#
# rule all:
#     input:
#         expand('results/JCASTLR_output/{sample}_JCASTLR{level}.fasta',sample=config['samples']['long'], level=["canonnical", "level1", "level2"])
#         expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
#         expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short'])


# load rules

include:
    "pt1/pt1.smk"
include:
    "pt2/pt2.smk"

# target rules
rule all:
    input:
        expand("results/qc/{qc_sample}.html",qc_sample=config['qc_list']),
        expand("results/qc/{qc_sample}_fastqc.zip",qc_sample=config['qc_list']),
        "results/qc/multiqc_report.html",
        expand("results/aligned/{sample}.{ext}",sample=config["samples"]["long"],ext=["bam", "bed"]),
        # expand("results/{sample}.productivity.bed", sample=config["samples"]),
        expand("{STAR_index}",STAR_index=config['STAR_index']),
        expand("results/aligned/{short_sample}_F_R_Aligned.sortedByCoord.out.bam",short_sample=config["samples"]["short"]),
        expand("results/junctions/{short_sample}_junctions.bed",short_sample=config['samples']['short']),
        expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
        expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
        expand("results/corrected/{sample}_all_corrected.bed", sample=config['samples']['long']),
        expand("results/corrected/{sample}_all_inconsistent.bed", sample=config['samples']['long']),
        expand("results/isoforms/{sample}.isoforms.{ext}", sample=config['samples']['long'], ext=["bed", "fa", "gtf"]),
        expand('results/JCASTLR_output/{sample}_JCASTLR{level}.fasta',sample=config['samples']['long'],level=["canonical", "level1", "level2"])
