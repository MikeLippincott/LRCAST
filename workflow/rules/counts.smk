
rule count_matrix_generation:
    input:
        bams = expand("results/aligned/long/{sample}.bam",sample=config['samples']['long']),
        gtf = "resources/genome/Homo_sapiens.GRCh38.107.gtf"
    output:
        out = "results/DGE/featurecounts.txt"
    log:
        "results/logs/DGE/featurecounts_log.txt"
    shell:
        "featureCounts -L -T 4 -a {input.gtf} -o {output.out} {input.bams} > {log}"

rule run_DESeq2:
    input:
        meta = "resources/meta.txt",
        counts = "results/DGE/featurecounts.txt"
    output:
        'DGE/normalized_counts.txt'
    shell:
        r"""
        scripts/DESeq2_DGE.R -i {input.counts} -m {input.meta} -o {output}
        """