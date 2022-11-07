from snakemake.utils import min_version
min_version("6.0")
configfile: "config/config.yaml"


include: "workflow/rules/qc.smk"
include: "workflow/rules/align_long.smk"
include: "workflow/rules/collapse_and_JCASTLR.smk"
include: "workflow/rules/counts.smk"
include: "workflow/rules/quantify.smk"

test = config["runmode"]

if test == "Long":
    include: "workflow/rules/correct_with_gtf.smk"
    rule all:
        input:
            # expand("results/qc/{qc_sample}.html",qc_sample=config['qc_list']),
            # expand("results/qc/{qc_sample}_fastqc.zip",qc_sample=config['qc_list']),\
            "results/qc/multiqc_report.html",
            expand("results/aligned/long/aligned.{ext}",sample=config["samples"]["long"],ext=["bam", "bed"]),
            "results/corrected/_all_corrected.bed",
            "results/corrected/_all_inconsistent.bed",

            expand("results/isoforms/collapse.isoforms.{ext}",ext=["bed", "fa","gtf"]),
            # expand('results/JCASTLR_output/{sample}_JCASTLR_{level}.fasta', \
            #     sample=config['samples']['long'],level=["Level1", "Level2", "Level3","Level4","Level5"]),
            # "results/DGE/featurecounts.txt",
            # 'results/DGE/normalized_counts.txt',
            # "results/qc/qualimap/qualimapReport.html",
            # "results/qc/qualimap/rnaseq_qc_results.txt"
            'results/DGE/counts_matrix.counts.tsv',
            'results/DGE/diffsplice.alt3.events.quant.tsv',
            'results/DGE/diffsplice.alt5.events.quant.tsv',
            'results/DGE/diffsplice.es.events.quant.tsv',
            'results/DGE/diffsplice.ir.events.quant.tsv',
            directory('results/DGE/diff_exp')

elif test == "Short":
    include: "workflow/rules/align_short.smk"
    include: "workflow/rules/correct_with_short.smk"
    include: "workflow/rules/counts.smk"
    rule all:
        input:
            expand("results/qc/{qc_sample}.html",qc_sample=config['qc_list']),
            expand("results/qc/{qc_sample}_fastqc.zip",qc_sample=config['qc_list']),\
            "results/qc/multiqc_report.html",
            expand("results/aligned/long/{sample}.{ext}",sample=config["samples"]["long"],ext=["bam", "bed"]),
            expand("{STAR_index}",STAR_index=config['STAR_index']),
            expand("results/aligned/short/{short_sample}_F_R_Aligned.sortedByCoord.out.bam",\
                short_sample=config["samples"]["short"]),
            expand("results/junctions/{short_sample}_junctions.bed",short_sample=config['samples']['short']),
            expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",zip,\
                sample=config['samples']['long'],short_sample=config['samples']['short']),
            expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",zip,\
                sample=config['samples']['long'],short_sample=config['samples']['short']),
            expand("results/corrected/{sample}_all_corrected.bed",sample=config['samples']['long']),
            expand("results/corrected/{sample}_all_inconsistent.bed",sample=config['samples']['long']),
            expand("results/isoforms/{sample}.isoforms.{ext}",\
                sample=config['samples']['long'],ext=["bed", "fa","gtf"]),
            expand('results/JCASTLR_output/{sample}_JCASTLR_{level}.fasta',\
                sample=config['samples']['long'],level=["Level1", "Level2", "Level3","Level4","Level5"]),
            "results/DGE/featurecounts.txt",
            'results/DGE/normalized_counts.txt',
            "results/qc/qualimap/qualimapReport.html",
            "results/qc/qualimap/rnaseq_qc_results.txt"
else:
    print("Error, please allow a proper run_type")
