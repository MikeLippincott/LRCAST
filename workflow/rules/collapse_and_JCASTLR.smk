rule run_flair_collapse:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = expand("resources/reads/{sample}.fastq",sample=config['samples']['long']),
        q = "results/corrected/_all_corrected.bed"
    output:
        bed = "results/isoforms/collapse.isoforms.bed",
        fa = "results/isoforms/collapse.isoforms.fa",
        gtf = "results/isoforms/collapse.isoforms.gtf",
    params:
        oname="results/isoforms/collapse"
    log:
        "results/logs/collapse/collapse.log"
    shell:
        """
        flair collapse -g {input.fa} -r {input.fq} -q {input.q} -f {input.gtf} -o {params.oname} 1>> {log} 2>>{log}
        """


