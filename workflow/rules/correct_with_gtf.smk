rule run_flair_correct_long:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        q = "results/aligned/long/aligned.bed",
    output:
        b1= "results/corrected/_all_corrected.bed",
        b2="results/corrected/_all_inconsistent.bed"
    log:
        "results/logs/corrected/corrected.log"
    params:
        oname = "results/corrected/"
    shell:
        """
        flair correct -q {input.q} -g {input.fa} -f {input.gtf} -o {params.oname} 1>> {log} 2>>{log}
        """
