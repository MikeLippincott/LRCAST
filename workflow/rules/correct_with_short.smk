rule run_flair_correct_short:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = "resources/reads/{sample}.fastq",
        q = "results/aligned/long/{sample}.bed",
        junc = "results/junctions/{short_sample}_junctions.bed"
    output:
        b1 ="results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",
        b2 = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",
    priority: 7
    params:
        oname = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}",
        oname1 = "results/corrected/{sample}_all_corrected.bed",
        oname2 = "results/corrected/{sample}_all_inconsistent.bed"
    log:
        "results/logs/corrected/{sample}_{short_sample}.log"
    shell:
        """
        flair correct -q {input.q} -g {input.fa} -f {input.gtf} -j {input.junc} -o {params.oname} && \
        sleep 5 && \
        cp {output.b1} {params.oname1} && \
        cp {output.b2} {params.oname2} && \
        sleep 15 1>> {log} 2>>{log}
        """



rule output_of_beds:
    output:
        expand("results/corrected/{sample}_all_corrected.bed",sample=config['samples']['long']),
        expand("results/corrected/{sample}_all_inconsistent.bed",sample=config['samples']['long'])