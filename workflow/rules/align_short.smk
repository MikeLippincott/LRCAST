rule star_index:
    input:
        fa = expand("{fa}",fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
    output:
        index = directory(f'{config["STAR_index"]}')
    params:
        sjdbOverhang = "149"
    log:
        "results/logs/star/index.log"
    shell:
        "STAR --runThreadN 4 --runMode genomeGenerate \
        --genomeDir {output.index} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdbOverhang} 1>> {log} 2>>{log}"

rule star_map:
    input:
        fq_F = "resources/reads/{short_sample}_1.fastq",
        fq_R = "resources/reads/{short_sample}_2.fastq",
        index = expand("{STAR_index}", STAR_index=config["STAR_index"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
    output:
        outfile = "results/aligned/short/{short_sample}_F_R_Aligned.sortedByCoord.out.bam"
    log:
        "results/logs/star/{short_sample}.log"
    params:
        outSAMtype = "BAM SortedByCoordinate",
        file_prefix = "results/aligned/short/{short_sample}_F_R_",
        overhang = "149"
    shell:
        "STAR --runThreadN 4 --genomeDir {input.index} \
        --sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} \
        --readFilesIn {input.fq_F} {input.fq_R} \
        --outSAMtype {params.outSAMtype} --outFileNamePrefix {params.file_prefix} 1>> {log} 2>>{log}"

rule run_flair_junctions:
    input:
        bam = "results/aligned/short/{short_sample}_F_R_Aligned.sortedByCoord.out.bam"
    output:
        junction = "results/junctions/{short_sample}_junctions.bed"
    log:
        "results/logs/junctions/{short_sample}.log"
    params:
        outname = "results/junctions/{short_sample}"
    shell:
        "junctions_from_sam -s {input.bam} -n {params.outname} 1>> {log} 2>>{log}"
