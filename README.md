# Long Read Centric Alternative Splicing Translator
## a.k.a
## mRNA centric Non-Junction Alternative Splicing Translator in silico
## McNJASTI


## Long
```
snakemake -j4 --ri --latency-wait 15
```

```
cd LRCAST
python JCAST -g /path/to/flair/gtf -f /path/to/flair/fa -p prefix/name -r cutoff_value -a altORF_Mode -o output/path
```

This workflow utilizes FLAIR and pieces from JCAST, and RIANA


