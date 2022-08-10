rule fastqc:
    input:
        "resources/reads/{qc_sample}.fastq"
    output:
        html="results/qc/{qc_sample}.html",
        zip="results/qc/{qc_sample}_fastqc.zip"
    # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "results/logs/qc/fastqc/{qc_sample}.log"
    threads: 4
    wrapper:
        "v1.3.2/bio/fastqc"

rule multiqc:
    input:
        expand("results/qc/{qc_sample}_fastqc.zip", qc_sample=config['qc_list'])
    output:
        "results/qc/multiqc_report.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "results/logs/qc/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"
