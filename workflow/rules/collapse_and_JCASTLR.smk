rule run_flair_collapse:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = "resources/reads/{sample}.fastq",
        q = "results/corrected/{sample}_all_corrected.bed"
    output:
        bed = "results/isoforms/{sample}.isoforms.bed",
        fa = "results/isoforms/{sample}.isoforms.fa",
        gtf = "results/isoforms/{sample}.isoforms.gtf",
    params:
        oname="results/isoforms/{sample}"
    log:
        "results/logs/collapse/{sample}.log"
    shell:
        """
        flair collapse -g {input.fa} -r {input.fq} -q {input.q} -f {input.gtf} -o {params.oname} 1>> {log} 2>>{log}
        """

rule run_JCAST_LR:
    input:
        gtf = "results/isoforms/{sample}.isoforms.gtf",
        fa = "results/isoforms/{sample}.isoforms.fa",
        bed = "results/isoforms/{sample}.isoforms.bed"
    output:
        'results/JCASTLR_output/{sample}_JCASTLR_{level}.fasta'
    params:
        prefix = "{sample}_",
        JCASTLR_path = "JCASTLR",
        out = "results/JCASTLR_output"
    log:
        "results/logs/JCASTLR/{sample}_{level}.log"
    shell:
        "python {params.JCASTLR_path} -g {input.gtf} -b {input.bed} -f {input.fa} -p {params.prefix} -o {params.out} 1>> {log} 2>>{log}"



