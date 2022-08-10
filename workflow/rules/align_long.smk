rule run_flair_align:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        fq = "resources/reads/{sample}.fastq"
    output:
        out = "results/aligned/long/{sample}.bam",
        out1 = "results/aligned/long/{sample}.bed"
    params:
        oname = "results/aligned/long/{sample}"
    log:
        "results/logs/flair_align/{sample}.log"
    shell:
        "flair align -g {input.fa} -r {input.fq} -o {params.oname} 1>> {log} 2>>{log}"
