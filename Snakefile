configfile: "resources/config.yaml"

rule all:
    input:
        # expand("results/qc/{sample}.html", sample=config["samples"]),
        # expand("results/qc/{sample}_fastqc.zip", sample=config["samples"]),
        expand("results/aligned/{sample}.{ext}", sample=config["samples"], ext=["sam","bam","bed"]),
        # "results/qc/multiqc_report.html",
        expand("results/{sample}_all_corrected.bed", sample=config["samples"]),
        expand("results/{sample}_all_inconsistent.bed", sample=config["samples"]),
        expand("results/{sample}.isoforms.{ext}", sample=config["samples"], ext=["bed","fa","gtf"])
        # expand("results/{sample}.productivity.bed", sample=config["samples"])


# rule fastqc:
#     input:
#         "resources/reads/{sample}.fastq.gz"
#     output:
#         html="results/qc/{sample}.html",
#         zip="results/qc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params: "--quiet"
#     log:
#         "results/logs/fastqc/{sample}.log"
#     threads: 4
#     wrapper:
#         "v1.3.2/bio/fastqc"
#
# rule multiqc:
#     input:
#         expand("results/qc/{sample}_fastqc.zip", sample=config["samples"])
#     output:
#         "results/qc/multiqc_report.html"
#     params:
#         ""  # Optional: extra parameters for multiqc.
#     log:
#         "results/logs/multiqc.log"
#     wrapper:
#         "v1.3.2/bio/multiqc"

# rule run_minimap:
#     input:
#         fa = "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#         gtf = "resources/genome/Homo_sapiens.GRCh38.106.gtf",
#         fq = expand("resources/reads/{sample}.fastq.gz", sample=config["samples"])
#     output:
#         sam = "results/aligned/{sample}.bam"
#     log:
#         "results/logs/{sample}.log"
#     shell:
#         "minimap2 -ax splice:hq -uf {input.fa} {input.fq} > {output.sam}"



rule run_flair_align:
    input:
        fa = "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf = "resources/genome/Homo_sapiens.GRCh38.106.gtf",
        fq = "resources/reads/{sample}.fastq.gz"
    output:
        out = "results/aligned/{sample}.{ext}"
    params:
        oname = "results/aligned/{sample}"
    # log:
        # "results/logs/{sample}.log"
    shell:
        "python flair/flair.py align -g {input.fa} -r {input.fq} -v1.3 -o {params.oname}"

rule run_flair_correct:
    input:
        fa = "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf = "resources/genome/Homo_sapiens.GRCh38.106.gtf",
        fq = "resources/reads/{sample}.fastq.gz",
        q = "results/aligned/{sample}.bed"
    output:
        b1 = "results/{sample}_all_corrected.bed",
        b2 = "results/{sample}_all_inconsistent.bed"
    params:
        oname="results/{sample}"
    # log:
    #     "results/logs/{sample}.log"
    shell:
        "python flair/flair.py correct -q {input.q} -g {input.fa} -f {input.gtf} -o {params.oname}"

rule run_flair_collapse:
    input:
        fa = "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf = "resources/genome/Homo_sapiens.GRCh38.106.gtf",
        fq = "resources/reads/{sample}.fastq.gz",
        q = "results/{sample}_all_corrected.bed"
    output:
        b1 = "results/{sample}.isoforms.bed",
        b2 = "results/{sample}.isoforms.fa",
        b3= "results/{sample}.isoforms.gtf"
    params:
        oname="results/{sample}",
        path="resources/reads/",
        samp="{sample}",
        fq = "resources/reads/{sample}.fastq"
    # log:
    #     "results/logs/{sample}.log"
    shell:
        """
        cp {input.fq} {params.path}/{params.samp}_1.fastq.gz
        gunzip {params.path}/{params.samp}_1.fastq.gz
        mv {params.path}/{params.samp}_1.fastq {params.path}/{params.samp}.fastq
        python flair/flair.py collapse -g {input.fa} -r {params.fq} -q {input.q} -f {input.gtf} -o {params.oname}
        """

# rule run_flair_predict:
#     input:
#         fa = "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#         gtf = "resources/genome/Homo_sapiens.GRCh38.106.gtf",
#         fq = "resources/reads/{sample}.fastq.gz",
#         q = "results/{sample}.isoforms.bed"
#     output:
#         b1 = "results/{sample}.productivity.bed"
#     # log:
#     #     "results/logs/{sample}.log"
#     shell:
#         "python flair/bin/predictProductivity.py -i {input.q} -g {input.gtf} -f {input.fa} --longestORF > {output.b1}"


