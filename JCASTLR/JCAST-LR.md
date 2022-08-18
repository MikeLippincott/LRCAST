

### Running JCAST-LR
```
python JCASTLR path/to/gtf/file.gtf -b path/to/bed/file.bed -f path/to/fasta/file.fa -o output/directory
```


```
>Canonnical|ENSG00000169972|ENST00000379031|J3KTG4 Q8N0Z8 |PUSL1|+|Chromosome 1 |  pseudouridine synthase like 1
```


JCAST-LR FASTA Custom Header Format 

1. Knowledgebase name, from canonical SwissProt protein entry (sp)
2. UniProt accession, from canonical SwissProt protein entry (Q91VW5)
3. UniProt name, from canonical SwissProt protein entry (GOGA4_MOUSE)
4. Annotated gene name (ENSMUSG00000038708)
5. Ensembl transcript id
6. Gene Symbol
7. Chromosome (chr9)
8. Isoform strand
9. Level (L1)


### JCAST-LR Levels:
### Level 1:
   Identified Gene ID and Transcript ID
### Level 2:
   Identified Gene ID Only
### Level 3:
   Orphaned Read