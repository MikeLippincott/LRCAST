
# rule count_matrix_generation:
#     input:
#         bams = 'results/aligned/long/aligned.bam',
#         gtf = "resources/genome/Homo_sapiens.GRCh38.107.gtf"
#     output:
#         out = "results/DGE/featurecounts.txt"
#     log:
#         "results/logs/DGE/featurecounts_log.txt"
#     shell:
#         "featureCounts -L -T 4 -a {input.gtf} -o {output.out} {input.bams} > {log}"

rule run_DESeq2:
    input:
        meta = "resources/experiment/meta.csv",
        counts = "results/DGE/diff_exp/filtered_gene_counts_ds2.tsv"
    output:
        'results/DGE/normalized_counts.txt'
    shell:
        """
        Rscript scripts/DESeq2_DGE.R -i {input.counts} -m {input.meta} -o {output}
        """