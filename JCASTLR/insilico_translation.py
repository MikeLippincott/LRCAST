# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Server # for retrieval of Uniprot IDs
import longread_functions as lrf


"""
Run this script with input files from the output of the FLAIR pipeline
script will output theorectical translated peptides from JCASTLR RNASeq data
"""



# Imports Bed File as object
class Bed:
    def __init__(self, bed_loc):
        self.bed = pd.read_csv(bed_loc, sep='\t', comment='t', header=None)
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                  'blockCount', 'blockSizes', 'blockStarts']
        self.bed.columns = header[:len(self.bed.columns)]
        # self.name = self.bed[self.bed['name'].str.contains('ENSG00000171163')]['name'].tolist()
        # self.score = self.bed[self.bed['name'].str.contains('ENSG00000171163')]['score'].tolist()
        # self.start = self.bed[self.bed['name'].str.contains('ENSG00000171163')]['chromStart'].tolist()
        # self.end = self.bed[self.bed['name'].str.contains('ENSG00000171163')]['chromEnd'].tolist()

# Imports GTF file as object
class Gtf:
    def __init__(self, gtf_loc):
        self.gtf = read_gtf(gtf_loc)
    def return_func(self):
        return self.gtf

# Loads info from biomart server for later queries
class LoadMart:
    def __init__(self):
        # Set Sever
        server = Server(host='http://www.ensembl.org')
        # retreive dataset used to query later
        self.dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])

# Sequence object
class Sequences(object):
    def __init__(self,
                 annotation_df: Gtf,
                 mart: LoadMart,
                 record):
        self.mart = mart.dataset
        self.gtf = annotation_df.return_func()
        self.transcript = record.id.split('_')[0]
        self.id = record.id.split('_')[1]
        self.seq = record.seq
    # Subset the GTF file for transcipt of interst
    def subset_gtf(self):
        self.gtf0 = self.gtf.query(f'gene_id == "{self.id}"').query('feature == "transcript"')
        self.strand = self.gtf0['strand'].unique()
        if len(self.strand) == 1:
            self.strand = self.strand[0]
        else:
            self.strand = 0
        self.frame = self.gtf0['frame'].unique()[0]
        if len(self.frame) == 1:
            self.frame = int(self.frame[0])
        else:
            self.frame == None
        # return self.id, self.strand, self.frame, self.transcript
    # Retrive meta data for transcript
    def get_meta(self):

        if 'ENST' in self.transcript:
            # print('Canonnical', self.transcript, self.id)

            enst = self.mart.query(attributes=["hgnc_symbol",
                                             'chromosome_name',
                                             'uniprot_gn_id',
                                             'entrezgene_description'],
                                 filters={'link_ensembl_transcript_stable_id': self.transcript})
            self.level = 'L1'
            self.gene_name = str(np.unique(enst['NCBI gene (formerly Entrezgene) description'])).strip("[']")
            self.gene_symbol = str(np.unique(enst['HGNC symbol'])).strip("[']")
            self.chromosome = str(np.unique(enst['Chromosome/scaffold name'])).strip("[']")

            if self.gene_name == float('nan'):
                self.gene_name = '---'
            if self.gene_symbol == float('nan'):
                self.gene_symbol = '---'
            if self.chromosome == float('nan'):
                self.chromosome = '---'

            self.uniprot = ''
            for i in np.unique(enst['UniProtKB Gene Name ID']):
                if str(type(i)) == "<class 'numpy.float64'>":
                    if np.isnan(i):
                        i = '---'
                i = str(i.strip("'"))
                self.uniprot += i
                self.uniprot += ' '
            if self.uniprot == float('nan'):
                self.uniprot = '---'
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        elif 'ENSG' in self.id:
            # print('Level 1', self.transcript, self.id)
            ensg = self.mart.query(attributes=["hgnc_symbol",
                                             'chromosome_name',
                                             'uniprot_gn_id',
                                             'entrezgene_description'], filters={'link_ensembl_gene_id': self.id})
            self.level = 'L2'
            self.gene_name = str(np.unique(ensg['NCBI gene (formerly Entrezgene) description'])).strip("[']")
            self.gene_symbol = str(np.unique(ensg['HGNC symbol'])).strip("[']")
            self.chromosome = str(np.unique(ensg['Chromosome/scaffold name'])).strip("[']")

            if self.gene_name == float('nan'):
                self.gene_name = '---'
            if self.gene_symbol == float('nan'):
                self.gene_symbol = '---'
            if self.chromosome == float('nan'):
                self.chromosome = '---'

            self.uniprot = ''
            for i in np.unique(ensg['UniProtKB Gene Name ID']):
                if str(type(i)) == "<class 'numpy.float64'>":
                    if np.isnan(i):
                        i = '---'
                i = str(i.strip("'"))
                self.uniprot += i
                self.uniprot += ' '
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        else:
            # print('Level 2', self.transcript,self.id)
            self.level = 'L3'
            self.gene_name = '---'
            self.gene_symbol = '---'
            self.uniprot = '---'
            self.chromosome = '---'
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome
    # make a header for outputed Fasta file in Peptide class
    # def make_header(self):
    #     self.header = "{0}|{1}|{2}|{3}|{4}|{5}|{6} | ".format(
    #         self.level,
    #         self.id,
    #         self.transcript,
    #         self.uniprot,
    #         self.gene_symbol,
    #         self.strand,
    #         f'Chromosome {self.chromosome}',)
    def make_header(self):
        self.header = "{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|".format(
            "sp",
            self.uniprot,
            self.id,
            self.transcript,
            self.gene_symbol,
            self.strand,
            f'Chromosome {self.chromosome}',
            self.level,)
        # return self.header

# Peptide Seq object class
class Peptide(object):
    def __init__(self,
                 sequence: Sequences):
        sequence.subset_gtf()
        sequence.get_meta()
        sequence.make_header()
        self.id = sequence.id
        self.seq = str(sequence.seq)
        self.strand = sequence.strand
        self.frame = sequence.frame
        self.header = sequence.header
        self.gene_name = sequence.gene_name

    # do a multiphase translation using the translate function to find ORFs
    def multi_phase_translate(self):
        """
        :return: longest peptide
        """
        d = {}
        for j in [0, 1, 2]:
            a = lrf.translate(self.seq, j)
            d[len(a)] = a
            self.prot = d[max(d)]
        # return self.prot
    # returns longest translated aa to a BioSeq Record object
    def str_to_seqrec(self):
        self.rec = SeqRecord(Seq(self.prot),f'{self.header}',description=self.gene_name)
        return self.rec

