rule run_flair_align:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        fq = expand("resources/reads/{sample}.fastq",sample=config['samples']['long'])
    output:
        out = "results/aligned/long/aligned.bam",
        out1 = "results/aligned/long/aligned.bed"
    params:
        oname = "results/aligned/long/aligned"
    log:
        "results/logs/flair_align/aligned.log"
    shell:
        "flair align -g {input.fa} -r {input.fq} -o {params.oname} 1>> {log} 2>>{log}"
