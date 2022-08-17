## Initialization
### Load and pasrse BED + GTF
### Load servers that need to accessed
call from insilico_translaion.py
## Retrieving Meta Data
call from insilico_translaion.py  
retrieve the following:  
- Name
- symbol
- Chromosome
- uniprot id
- define level
- Biotype
- support level (tsl)
- Number of Exons

# Define Levels
### Level 1
ENSG, ENST, Protein Coding
### Level 2
ENSG, protein coding
### Level 3
ENSG, ENST, non-protein coding  
ENSG, non-protein coding  
### Level 4
No annotation found


## Translation
### Create Header for output Fasta
Logic:  
### If Level 1  
translate annotated transcript from TSS  
### If Level 2   
translate from cannonical TSS  
&  
translate multiphase  
### If Level 3  
translate multiphase  
### If Level 4   
translate multiphase  

Levels 3 and 4 set peptide length cut-off for micropeptide DB

## Final Piece
write peptides to fasta db  

Add Support for exon usage annotation
Add support for uORF annotaion
Add seperate function for riboseq data integration





