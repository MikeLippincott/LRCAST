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

rule run_qualimap:
    input:
        gtf = "resources/genome/Homo_sapiens.GRCh38.107.gtf",
        bams = expand("results/aligned/long/{sample}.bam",sample=config['samples']['long'])
    output:
        o1 = "results/qc/qualimap/qualimapReport.html",
        o2 = "results/qc/qualimap/rnaseq_qc_results.txt"
    params:
        outpath = "results/qc/qualimap/",
        mem = "8G"
    shell:
        "qualimap rnaseq -outdir {params.outpath} -a proportional \
        -bam {input.bams} -gtf {input.gtf} --java-mem-size={params.mem}"
