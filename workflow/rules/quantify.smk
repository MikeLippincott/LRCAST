rule run_flair_quantify:
    input:
        fa = expand("{fa}", fa=config["genome"]["fa"]),
        gtf = expand("{gtf}", gtf=config["genome"]["gtf"]),
        isos = 'results/isoforms/collapse.isoforms.fa',
        reads = "resources/experiment/experiment_info.tsv"
    output:
        b1= 'results/DGE/counts_matrix.counts.tsv'
    log:
        "results/logs/quantify.log"
    params:
        oname = "results/DGE/counts_matrix"
    shell:
        """
        flair quantify -r {input.reads} -i {input.isos} -o {params.oname} 1>> {log} 2>>{log}
        """

rule run_flair_diffsplice:
    input:
        q = "results/DGE/counts_matrix.counts.tsv",
        isos= 'results/isoforms/collapse.isoforms.bed',
    output:
        'results/DGE/diffsplice.alt3.events.quant.tsv',
        'results/DGE/diffsplice.alt5.events.quant.tsv',
        'results/DGE/diffsplice.es.events.quant.tsv',
        'results/DGE/diffsplice.ir.events.quant.tsv'
    params:
        oname = "results/DGE/diffsplice"
    shell:
        'flair diffSplice -i {input.isos} -q {input.q} -o {params.oname}'


rule run_flair_diff_exp:
    input:
        q="results/DGE/counts_matrix.counts.tsv"
    output:
        directory('results/DGE/diff_exp/')
    params:
        oname = "results/DGE/diff_exp"
    shell:
        'flair diffExp -q {input.q} -o {params.oname}'


