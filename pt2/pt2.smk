import os
# configfile: "resources/config.yaml"

# rule all:
#     input:
        # expand("results/corrected/{sample}_all_corrected.bed", sample=config['samples']['long']),
        # expand("results/corrected/{sample}_all_inconsistent.bed", sample=config['samples']['long']),
        # expand("results/logs/tmp/{short_sample}_{sample}.txt", sample=config['samples']['long'],short_sample=config['samples']['short']),
        # expand("results/isoforms/{sample}.isoforms.{ext}", sample=config['samples']['long'], ext=["bed", "fa", "gtf"]),
        # expand('results/JCASTLR_output/{sample}_JCASTLR{level}.fasta',sample=config['samples']['long'], level=["canonnical", "level1", "level2"])

rule run_flair_collapse:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = "resources/reads/{sample}.fastq",
        # q = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed"
        q = "results/corrected/{sample}_all_corrected.bed"
    output:
        bed = "results/isoforms/{sample}.isoforms.bed",
        fa = "results/isoforms/{sample}.isoforms.fa",
        gtf = "results/isoforms/{sample}.isoforms.gtf",
    params:
        oname="results/isoforms/{sample}"
    priority: 8
    log:
        "results/logs/{sample}.log"
    shell:
        """
        flair collapse -g {input.fa} -r {input.fq} -q {input.q} -f {input.gtf} -o {params.oname}
        """

rule run_JCAST_LR:
    input:
        # gtf = rules.run_flair_collapse.output.gtf,
        # fa = rules.run_flair_collapse.output.fa,
        # bed = rules.run_flair_collapse.output.bed
        gtf = "results/isoforms/{sample}.isoforms.gtf",
        fa = "results/isoforms/{sample}.isoforms.fa",
        bed = "results/isoforms/{sample}.isoforms.bed"
    output:
        'results/JCASTLR_output/{sample}_JCASTLR{level}.fasta'
    params:
        prefix = "{sample}_",
        JCASTLR_path = "JCASTLR",
        out = "results/JCASTLR_output"
    priority: 9
    shell:
        "python {params.JCASTLR_path} -g {input.gtf} -b {input.bed} -f {input.fa} -p {params.prefix} -o {params.out}"
