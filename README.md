# Long Read Centric Alternative Splicing Translator
## a.k.a
## mRNA centric Non-Junction Alternative Splicing Translator in silico
## McNJASTI


#### LRCAST uses output files from the flair pipeline (Brooks Lab) to generate a protein Database

### Run snakemake for flair pipeline processing 
```
snakemake --cores all --ri --latency-wait 15
```

### Run LRCAST
```
cd LRCAST
python JCAST -g /path/to/flair/gtf -f /path/to/flair/fa -p prefix/name -r cutoff_value<integer> -o output/path
```

##### Additional Args: 
```
-a #altORF generation mode 
```

This workflow utilizes FLAIR and pieces from JCAST, and RIANA


