rule run_flair_correct_long:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = "resources/reads/{sample}.fastq",
        q = "results/aligned/long/{sample}.bed",
    output:
        b1= "results/corrected/{sample}_all_corrected.bed",
        b2="results/corrected/{sample}_all_inconsistent.bed"
    log:
        "results/logs/corrected/{sample}.log"
    params:
        oname = "results/corrected/{sample}"
    shell:
        """
        flair correct -q {input.q} -g {input.fa} -f {input.gtf} -o {params.oname} 1>> {log} 2>>{log}
        """
