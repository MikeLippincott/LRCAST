configfile: "resources/config.yaml"


# rule all:
#     input:
#         expand("results/qc/{qc_sample}.html",qc_sample=config['qc_list']),
#         expand("results/qc/{qc_sample}_fastqc.zip",qc_sample=config['qc_list']),
#         "results/qc/multiqc_report.html",
#         expand("results/aligned/{sample}.{ext}", sample=config["samples"]["long"], ext=["bam","bed"]),
#         # expand("results/{sample}.productivity.bed", sample=config["samples"]),
#         expand("{STAR_index}", STAR_index=config['STAR_index']),
#         expand("results/aligned/{short_sample}_F_R_Aligned.sortedByCoord.out.bam", short_sample=config["samples"]["short"]),
#         expand("results/junctions/{short_sample}_junctions.bed", short_sample=config['samples']['short']),
#         expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed", zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
#         expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed", zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
#         # expand("results/corrected/{sample}_all_corrected.bed", sample=config['samples']['long']),
        # expand("results/corrected/{sample}_all_inconsistent.bed", sample=config['samples']['long']),
        # expand("results/logs/tmp/{short_sample}_{sample}.txt", sample=config['samples']['long'],short_sample=config['samples']['short']),
        # expand("results/isoforms/{sample}.isoforms.{ext}", sample=config['samples']['long'], ext=["bed", "fa", "gtf"]),
        # expand('results/JCASTLR_output/{sample}_JCASTLR{level}.fasta',sample=config['samples']['long'], level=["canonnical", "level1", "level2"])

rule fastqc:
    input:
        "resources/reads/{qc_sample}.fastq"
    output:
        html="results/qc/{qc_sample}.html",
        zip="results/qc/{qc_sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    priority: 1
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
    priority: 2
    log:
        "results/logs/qc/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"

rule run_flair_align:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        fq = "resources/reads/{sample}.fastq"
    output:
        out = "results/aligned/{sample}.{ext}"
    params:
        oname = "results/aligned/{sample}"
    priority: 3
    # log:
    #     "results/logs/flair_align/{sample}.log"
    shell:
        "flair align -g {input.fa} -r {input.fq} -o {params.oname}"

rule star_index:
    input:
        fa = expand("{fa}",fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
    output:
        index = directory(f'{config["STAR_index"]}')
    params:
        sjdbOverhang = "149"
    priority: 4
    shell:
        "STAR --runThreadN 4 --runMode genomeGenerate \
        --genomeDir {output.index} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdbOverhang}"

rule star_map:
    input:
        fq_F = "resources/reads/{short_sample}_1.fastq",
        fq_R = "resources/reads/{short_sample}_2.fastq",
        index = expand("{STAR_index}", STAR_index=config["STAR_index"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
    output:
        outfile = "results/aligned/{short_sample}_F_R_Aligned.sortedByCoord.out.bam"
    priority: 5
    params:
        outSAMtype = "BAM SortedByCoordinate",
        file_prefix = "results/aligned/{short_sample}_F_R_",
        overhang = "149"
    shell:
        "STAR --runThreadN 4 --genomeDir {input.index} \
        --sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} \
        --readFilesIn {input.fq_F} {input.fq_R} \
        --outSAMtype {params.outSAMtype} --outFileNamePrefix {params.file_prefix}"

rule run_flair_junctions:
    input:
        bam = "results/aligned/{short_sample}_F_R_Aligned.sortedByCoord.out.bam"
    output:
        junction = "results/junctions/{short_sample}_junctions.bed"
    priority: 6
    params:
        outname = "results/junctions/{short_sample}"
    shell:
        "junctions_from_sam -s {input.bam} -n {params.outname}"

rule run_flair_correct:
    input:
        # bam = "results/aligned/{short_sample}_F_R_Aligned.sortedByCoord.out.bam",
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        fq = "resources/reads/{sample}.fastq",
        q = "results/aligned/{sample}.bed",
        junc = "results/junctions/{short_sample}_junctions.bed"
        # junc = expand("results/junctions/{short_sample}_junctions.bed", short_sample=config["samples"]["short"])
        # junc = lambda wildcards: ["results/junctions/{short_sample}_junctions.bed".format(short_sample) for short_sample in config["samples"]["short"]]
        # junc = [x for x in config["samples"]["short"]]
    output:
        b1 ="results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",
        b2 = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",
        # junction = "results/junctions/{short_sample}_junctions.bed",
        # b3= "results/corrected/{sample}_all_corrected.bed",
        # b4="results/corrected/{sample}_all_inconsistent.bed"
    priority: 7
    params:
        oname = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}",
        # oname = "results/corrected/{sample}",
        # b1 ="results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",
        # b2 = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",
        oname1 = "results/corrected/{sample}_all_corrected.bed",
        oname2 = "results/corrected/{sample}_all_inconsistent.bed"
    # log:
    #     "results/logs/{sample}.log"
    shell:
        """
        flair correct -q {input.q} -g {input.fa} -f {input.gtf} -j {input.junc} -o {params.oname} \
        && \
        cp {output.b1} {params.oname1} \
        && \
        cp {output.b2} {params.oname2}
        """

rule output_of_beds:
    input:
        "results/aligned/{sample}.bed"
    output:
        b1 = "results/corrected/{sample}_all_corrected.bed",
        b2 = "results/corrected/{sample}_all_inconsistent.bed"
#
# rule rename_files:
#     input:
#         b1 = rules.run_flair_correct.output.b1,
#         b2 = rules.run_flair_correct.output.b2
#         # b1 = expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
#         # b2 = expand("results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_inconsistent.bed",zip,sample=config['samples']['long'],short_sample=config['samples']['short']),
#     output:
#         b3 = "results/corrected/{sample}_all_corrected.bed",
#         b4 = "results/corrected/{sample}_all_inconsistent.bed",
#         # b5 = "results/logs/tmp/{sample}_{short_sample}_all_corrected.txt",
#         # b6 = "results/logs/tmp/{sample}_{short_sample}all_inconsistent.txt"
#     priority: 8
#     shell:
#         """
#         cp {input.b1} {output.b3} \
#         && \
#         cp {input.b2} {output.b4} \
#         """

# rule run_flair_collapse:
#     input:
#         fa = expand("{fa}", fa=config["genome"]["fa"]),
#         gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
#         fq = "resources/reads/{sample}.fastq",
#         # q = "results/corrected/flair_correct/{sample}_correctedby_{short_sample}_all_corrected.bed"
#         q= "results/corrected/flair_correct/{sample}_all_corrected.bed"
#     output:
#         "results/isoforms/{sample}.isoforms.bed",
#         "results/isoforms/{sample}.isoforms.fa",
#         "results/isoforms/{sample}.isoforms.gtf",
#     params:
#         oname="results/isoforms/{sample}"
#     priority: 8
#     log:
#         "results/logs/{sample}.log"
#     shell:
#         """
#         flair collapse -g {input.fa} -r {input.fq} -q {input.q} -f {input.gtf} -o {params.oname}
#         """
#
# # # ruleorder: run_flair_correct > rename_files
# rule run_JCAST_LR:
#     input:
#         gtf = "results/isoforms/{sample}.isoforms.gtf",
#         fa = "results/isoforms/{sample}.isoforms.fa",
#         bed = "results/isoforms/{sample}.isoforms.bed"
#     output:
#         'results/JCASTLR_output/{sample}_JCASTLR{level}.fasta'
#     params:
#         prefix = "{sample}_",
#         out = "results/JCASTLR_output"
#     priority: 9
#     shell:
#         "python JCASTLR -g {input.gtf} -b {input.bed} -f {input.fa} -p {params.prefix} -o {params.out}"


