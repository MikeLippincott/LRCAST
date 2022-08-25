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


rule count_matrix_generation:
    input:
        bams = expand("results/aligned/long/{sample}.bam",sample=config['samples']['long']),
        gtf = "resources/genome/Homo_sapiens.GRCh38.107.gtf"
    output:
        out = "results/DGE/featurecounts.txt"
    log:
        "results/logs/DGE/featurecounts.txt"
    shell:
        "featureCounts -L -T 4 -a {input.gtf} -o {output.out} {input.bams} > {log}"
