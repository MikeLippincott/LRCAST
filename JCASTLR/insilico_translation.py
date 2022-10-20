# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Server # for retrieval of Uniprot IDs
import constants
import longread_functions as lrf
# from JCASTLR import longread_functions as lrf
import gget
import sqlite3 as sq
import re
import os
import logging
import requests as rq
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from io import StringIO
from time import sleep
import model
import logging

"""
Run this module with input files from the output of the FLAIR pipeline
script will output theorectical translated peptides from JCASTLR RNASeq data
"""

# Imports GTF file as object
class Gtf:
    def __init__(self, gtf_loc, counts_matrix):
        print("LOADING GTF FILES INTO PANDAS DATAFRAMES.")
        # FLAIR output gtf
        self.gtf = read_gtf(gtf_loc)
        # Ensembl gtf
        self.gtf_file = read_gtf("resources/genome/Homo_sapiens.GRCh38.107.gtf")
        self.counts = pd.read_table(counts_matrix)

    def read_cutoff(self, write_dir):
        exp = pd.read_table('resources/experiment/experiment_info.tsv', header=None)
        lst = []
        for i, j, k in zip(exp[0], exp[1], exp[2]):
            col = f'{i}_{j}_{k}'
            lst.append(col)
        count = 0
        for i in lst:
            count += self.counts[f'{i}']
        self.counts['counts'] = count
        df1 = self.counts[(self.counts.counts > 1)]
        write_dir = os.getcwd()
        ln_sjc, best_mix_model, self.min_count = model.general_mixture_model(sum_sjc_array=df1['counts'].astype(int))
        model.plot_general_mixture_model(
            ln_sjc,
            best_mix_model,
            self.min_count,
            write_dir=write_dir,
            filename='model',
        )
        df3 = self.counts[(self.counts.counts > self.min_count)]
        self.gtf['ids'] = self.gtf["transcript_id"] + "_" + self.gtf["gene_id"]
        self.gtf = pd.merge(self.gtf, df3, on=['ids'], how='right')
        gene_ids = np.unique(self.gtf['gene_id'])
        print(f'{len(gene_ids)} Filtered transcripts to Process.')
        bools = self.gtf_file.gene_id.isin(gene_ids)
        self.gtf_file = self.gtf_file[bools]



# Sequence object
class Sequences(object):
    def __init__(self,
                 annotation_df: Gtf,
                 record):
        # FLAIR GTF
        self.gtf = annotation_df.gtf
        # Genomic GTF
        self.gtf_file = annotation_df.gtf_file
        self.counts_df = annotation_df.counts
        self.rid = record.id
        self.tid = self.rid.split('_')[0]
        self.gid = self.rid.split('_')[1]
        self.seq = record.seq

    # Subset the GTF file for transcipt of interst
    def subset_gtf(self):
        """
        queries the ensembl gtf for all rows belonging to gene and retrieves strand/frame (phase) information
        :return: no return
        """
        self.gtf0 = self.gtf.query(f'gene_id == "{self.gid}"').query('feature == "transcript"')
        self.strand = self.gtf0['strand'].unique()
        if len(self.strand) == 1:
            self.strand = self.strand[0]
        else:
            self.strand = 0
        self.frame = self.gtf0['frame'].unique()
        if len(self.frame) == 1:
            self.frame = int(self.frame[0])
        else:
            self.frame == None
        # return self.id, self.strand, self.frame, self.transcript

    def get_counts(self):
        # print(self.tid,self.gid)
        exp = pd.read_table('resources/experiment/experiment_info.tsv', header=None)
        df1 = self.counts_df.loc[self.counts_df['ids'] == f'{self.tid}_{self.gid}']
        lst = []
        for i, j, k in zip(exp[0], exp[1], exp[2]):
            col = f'{i}_{j}_{k}'
            lst.append(col)
        t = 0
        for i in lst:
            if len(df1[f'{i}']) == 0:
                self.counts = 0
            else:
                t += int(df1[f'{i}'])
        self.counts = t
        # print(self.counts)


    # retrieve transcript or gene meta data
    def get_meta(self):

        """
        Function defines each record's level then retrieves transcript or gene meta data depending on which ensembl id
        is available (ensembl transcript id has priority over ensembl gene id)
        Metadata to retieve:
        level
        biotype
        gene symbol
        gene name
        chromosome
        transcript support level (tsl)
        uniprot id
        :return: No return
        """
        # print(self.rid)
        if 'ENST' in self.rid:
            enst = self.gtf_file.query(f'transcript_id == "{self.tid}"').query('feature == "transcript"')
            gtf0 = self.gtf_file.query(f'transcript_id == "{self.tid}"').query('transcript_biotype == "protein_coding" ')
            if len(gtf0) > 0:
                self.level = "L1"
                self.biotype = "protein_coding"
            elif len(gtf0) == 0:
                self.level = "L3"
                biotype = ''
                tmp = np.unique(self.gtf_file.query(f'transcript_id == "{self.tid}"')['transcript_biotype'])
                if len(tmp) > 0:
                    if len(biotype) > 0:
                        biotype += '__'
                    else:
                        for i in range(len(tmp)):
                            if tmp[i] == '':
                                biotype = biotype
                            else:
                                biotype += tmp[i]
                else:
                    biotype = '-'
                self.biotype = biotype
            else:
                print("Check Meta Data method in Sequences Class")

            sleep(0.1)
            if pd.isnull(enst['gene_name'].to_list()):
                self.gene_symbol = '-'
            elif enst['gene_name'].to_list() == []:
                self.gene_symbol = '-'
            else:
                self.gene_symbol = enst['gene_name'].to_list()[0]

            self.chromosome = self.gtf0['seqname'].to_list()[0]
            if len(gtf0['transcript_support_level']) == 0:
                self.tsl = "-"
            else:
                self.tsl = gtf0['transcript_support_level'].to_list()[0]
            # self.uniprot = enst['uniprot_id'].to_list()[0]
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        elif 'ENSG' in self.rid:
            gtf0 = self.gtf_file.query(f'gene_id == "{self.gid}"').query('transcript_biotype == "protein_coding" ')
            if len(gtf0) > 0:
                self.level = "L2"
                self.biotype = "protein_coding"
            elif len(gtf0) == 0:
                biotype = ''
                tmp = np.unique(self.gtf_file.query(f'gene_id == "{self.gid}"')['transcript_biotype'])
                if len(tmp) > 0:
                    if len(biotype) > 0:
                        biotype += '__'
                    else:
                        for i in range(len(tmp)):
                            if tmp[i] == '':
                                biotype = biotype
                            else:
                                biotype += tmp[i]
                else:
                    biotype = '-'
                self.biotype = biotype


                self.level = "L4"
            else:
                print("Check Meta Data method in Sequences Class")

            # print(self.gid)
            ensg = self.gtf_file.query(f'gene_id == "{self.gid}"').query('feature == "gene"')
            sleep(0.1)

            # print(self.gene_name)
            if pd.isnull(ensg['gene_name'].to_list()):
                self.gene_symbol = '-'
            elif ensg['gene_name'].to_list() == []:
                self.gene_symbol = '-'
            else:
                self.gene_symbol = ensg['gene_name'].to_list()[0]
            self.chromosome = self.gtf0['seqname'].to_list()[0]
            # self.uniprot = ensg['uniprot_id'].to_list()[0]
            if len(gtf0['transcript_support_level']) == 0 :
                self.tsl = "-"
            else:
                self.tsl = gtf0['transcript_support_level'].to_list()[0]
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        elif not "ENSG" in self.rid:
            self.level = "L5"
            self.biotype = "L5"
            # self.gene_name = "L5"
            self.gene_symbol = "L5"
            self.chromosome = self.gtf0['seqname'].to_list()[0]
            # self.uniprot = "L5"
            self.tsl = "L5"

        else:
            print('LRCAST: RNA Fasta Header Read Error!')

        # if pd.isnull(self.gene_name):
        #     self.gene_name = '-'

    def annotated_trancript_trim(self):
        """
        :return: trimmed sequence for insilico translation
        """
        gtf0 = self.gtf_file.query(f'transcript_id == "{self.tid}"').query('transcript_biotype == "protein_coding" ')
        TSS = gtf0.query('feature == "start_codon" ')
        if TSS['start'].to_list() == []:
            gtf0 = self.gtf_file.query(f'gene_id == "{self.gid}"').query('transcript_biotype == "protein_coding" ')
            TSS = gtf0.query('feature == "start_codon" ')
        transcript_df = gtf0.query('feature == "transcript" ')
        exons = gtf0.query('feature == "exon" ')
        # exons = exons.sort_values(by="exon_number")
        jcast = self.gtf.query(f'transcript_id == "{self.tid}"').query('feature == "transcript" ')
        x1 = jcast['start'].to_list()[0]
        x2 = jcast['end'].to_list()[0]
        # print(transcript_id)
        # print(TSS['strand'].to_list())

        strand = jcast['strand'].to_list()[0]
        phase = int(jcast['frame'].to_list()[0])
        if len(TSS) == 0:
            print("No Annotated Start Codon")
            self.a_seq = ''
            return self.a_seq
        elif len(TSS) > 0:
            y1 = TSS['start'].to_list()[0]
            y2 = TSS['end'].to_list()[0]
        if strand == '+':
            SS = y1 - x1
        elif strand == '-':
            SS = x2 - y2
        else:
            print("Start Site Extraction Error")
        # print(f'TSS Coordinate: {SS}')
        # print(f'The Length of This Transcript is, {-(x1 - x2)}')
        start = exons['start'].to_list()
        end = exons['end'].to_list()
        i = 0
        e = 0
        while i <= len(end) - 1:
            if start[i] < y1 < end[i]:
                # print(start[i],", ",y1,", ", end[i])
                e = i + 1
                # print(f'start Codon for {transcript_id} is in Exon {e}')
                i += 1
            else:
                i += 1
        if strand == "+":
            start = start[:e]
            end = end[:e]
            i = 1
            j = 0
            t_len = 0
            while i <= e - 1:
                t_len += (start[i] - end[j]) - 1
                j += 1
                i += 1
            adjustedSS = (SS - t_len)
            # print(adjustedSS)
        elif strand == "-":
            start = start[:e]
            end = end[:e]
            i = 0
            j = 1
            t_len = 0
            while j <= e - 1:
                t_len += (start[i] - end[j]) - 1
                j += 1
                i += 1
            adjustedSS = (SS - (t_len))
        self.a_seq = (str(self.seq[adjustedSS:]))
        return self.a_seq


# Peptide Seq object class

class Peptide(object):
    def __init__(self,
                 sequence: Sequences):
        # sequence.subset_gtf()
        # sequence.get_meta()
        # sequence.annotated_trancript_trim()
        self.s = sequence
        self.sequence = sequence
        self.seq = str(sequence.seq)
        self.strand = sequence.strand
        self.frame = sequence.frame
        # self.gene_name = sequence.gene_name
        self.level = sequence.level
        self.a_seq = ''

    def annotated_translate(self):
        """
        used the trimmed sequence to translate from start codon
        :return: peptide from annotation coordinates
        """
        code = constants.genetic_code
        self.prot = ''
        # print(self.s.a_seq)
        if self.s.a_seq == '':
            self.prot = ''
        else:
            for i in range(0,len(self.s.a_seq) - 2, 3):
                if 'N' in self.sequence.a_seq[i:i+3]:
                    aa = 'U'
                else:
                    aa = code[self.s.a_seq[i:i + 3]]
                if aa == 'X':
                    return self.prot
                else:
                    self.prot += aa
            # print(self.prot)
        return self.prot

    # @staticmethod
    # def start_find_translate(nt, phase):
    #     """
    #     :nt: input nucleotide sequence
    #     :phase: integer 0-2 which phase to translate in
    #     :return: longest peptide
    #     """
    #     code = constants.genetic_code
    #     pep = ''
    #     for i in range(phase, len(nt) - 2, 3):
    #         if 'N' in nt[i:i+3]:
    #             aa = 'U'
    #         else:
    #             aa = code[nt[i:i + 3]]
    #         if len(pep) == 0:
    #             if aa == 'M':
    #                 pep += aa
    #                 # print(i)
    #         elif len(pep) > 0:
    #             if aa == 'X':
    #                 return pep
    #             else:
    #                 pep += aa
    #     # print(len(pep),", ",pep)
    #     return pep

    @staticmethod
    def start_find_translate(nt, phase):
        """
        :nt: input nucleotide sequence
        :phase: integer 0-2 which phase to translate in
        :return: longest peptide
        """
        code = constants.genetic_code
        pep = ''
        start = 0
        stop = 0
        for i in range(phase, len(nt) - 2, 3):
            if 'N' in nt[i:i + 3]:
                aa = 'U'
            else:
                aa = code[nt[i:i + 3]]
            if len(pep) == 0:
                if aa == 'M':
                    start = i
                    pep += aa
                    # print(i)
            elif len(pep) > 0:
                if aa == 'X':
                    stop = i
                    if start == 0:
                        start = 0
                    if stop == 0:
                        stop = 0
                    return pep, start, stop
                else:
                    pep += aa
        # print(len(pep),", ",pep)
        return pep, start, stop

    # do a multiphase translation using the translate function to find ORFs
    # def multi_phase_translate(self):
    #     """
    #     :return: longest peptide
    #     """
    #     d = {}
    #     for j in [0, 1, 2]:
    #         a = Peptide.start_find_translate(self.seq, j)
    #         d[len(a)] = a
    #         self.prot = d[max(d)]
    #     return self.prot

    # def multi_phase_translate(self):
    #     """
    #     :nt: input nucleotide sequence
    #     :return: longest peptide
    #     """
    #     t_len = 0
    #     d0 = {}
    #     d1 = {}
    #     d2 = {}
    #     fd = {}
    #     lst0 = []
    #     lst1 = []
    #     lst2 = []
    #     while len(self.seq) > t_len:
    #         t_len += 3
    #
    #         for j in [0, 1, 2]:
    #             # a = start_find_translate(seq, j)
    #             # d[a] = len(a)
    #             # prot = [max(d)]
    #
    #             a = Peptide.start_find_translate(self.seq[t_len:], j)
    #             # print(a)
    #             if j == 0:
    #                 lst0.append(a)
    #             elif j == 1:
    #                 lst1.append(a)
    #             elif j == 2:
    #                 lst2.append(a)
    #
    #     for k in lst0:
    #         d0[k] = len(k)
    #     for k in lst1:
    #         d1[k] = len(k)
    #     for k in lst2:
    #         d2[k] = len(k)
    #     for i, j in zip([d0, d1, d2], range((3))):
    #         fd[f'Phase{j}'] = i
    #
    #     self.prot = 'M'
    #     for i in fd:
    #         # print(fd[i])
    #         for j in fd[i]:
    #             if len(j) > len(self.prot):
    #                 self.prot = j
    #                 p = i
    #             else:
    #                 pass
    #
    #     return self.prot

    def multi_phase_translate(self):
        """
        :nt: input nucleotide sequence
        :return: longest peptide
        """
        t_len = 0
        d0 = {}
        d1 = {}
        d2 = {}
        fd = {}
        lst0 = []
        lst1 = []
        lst2 = []
        total_length = len(self.seq)
        start = 0
        stop = 0
        self.dict2 = {}
        while len(self.seq) > t_len:
            t_len += 3
            dict1 = {'phase': [], 'len': [], 'start': [], 'stop': []}
            for j in [0, 1, 2]:
                # a = start_find_translate(seq, j)
                # d[a] = len(a)
                # prot = [max(d)]

                a, start, stop = Peptide.start_find_translate(self.seq[t_len:], j)
                dict1 = {'phase': [], 'len': [], 'start': [], 'stop': []}
                # print(a)
                if j == 0:
                    dict1['phase'].append(f'Phase{j}')
                    dict1['len'].append(len(a))
                    dict1['start'].append(start + t_len)
                    dict1['stop'].append(stop + t_len)
                    # dict2[a] = dict1
                elif j == 1:
                    dict1['phase'].append(f'Phase{j}')
                    dict1['len'].append(len(a))

                    dict1['start'].append(start + t_len)
                    dict1['stop'].append(stop + t_len)
                    # dict2[a] = dict1

                elif j == 2:
                    dict1['phase'].append(f'Phase{j}')
                    dict1['len'].append(len(a))
                    dict1['start'].append(start + t_len)
                    dict1['stop'].append(stop + t_len)
                self.dict2[a] = dict1

        # for k in lst0:
        #     d0[k] = len(k)
        # for k in lst1:
        #     d1[k] = len(k)
        # for k in lst2:
        #     d2[k] = len(k)
        # for i,j in zip([d0,d1,d2],range((3))):
        #     fd[f'Phase{j}'] = i

        self.prot = 'M'
        for j in self.dict2:
            # print(fd[i])
            if len(j) > len(self.prot):
                self.prot = j
            else:
                pass

        # return k
        return self.prot

        # def get_canonical_aa_uniprot_local(self,
    #                                    ) -> SeqRecord:
    #     self.canonical_aa = ''
    #     old = self.level
    #     with open('resources/DB/reviewed_canonical.fasta') as f:
    #         for record in SeqIO.parse(f, 'fasta'):
    #             if record.seq == self.prot:
    #                 self.level = 'Canonical'
    #                 self.canonical_aa = record.id
    #     print(old,' ',self.level)
    #     print(self.canonical_aa)

    def make_header(self):
        self.header = "{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|".format(
            self.s.level,
            self.s.gene_symbol,
            self.s.gid,
            self.s.tid,
            self.s.strand,
            f'Chr{self.s.chromosome}',
            self.s.biotype,
            self.s.tsl,
            self.s.counts,)


    # returns longest translated aa to a BioSeq Record object
    def str_to_seqrec(self, prot_seq):
        self.rec = SeqRecord(Seq(prot_seq),f'{self.header}',description=self.s.gene_symbol)
        return self.rec



class Post_hoc_reassignment():
    def __init__(self,
                 sequence:Sequences,
                 peptide: Peptide):
        self.s = sequence
        self.p = peptide
        self.level = self.s.level
        self.header = self.p.header
        self.id = '-'


    def get_canonical_aa_uniprot_local(self,
                                       DB,
                                       ) -> SeqRecord:
        self.canonical_aa = ''
        self.old = self.level
        with open(DB) as f:
            for record in SeqIO.parse(f, 'fasta'):
                if record.seq == self.p.prot:
                    self.level = 'Canonical'
                    self.id = record.id
        # print(self.old, ' ', self.level)

    def get_aa_uniprot_local(self,
                             DB,
                             ) -> SeqRecord:
        self.canonical_aa = ''
        self.old = self.level
        with open(DB) as f:
            for record in SeqIO.parse(f, 'fasta'):
                if record.seq == self.p.prot:
                    self.id = record.id
        # print(self.old, ' ', self.level)
        # print(self.canonical_aa)


    def make_header(self):
        if self.id != '-':
            self.header = self.id + "|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|".format(
                self.s.level,
                self.s.gene_symbol,
                self.s.gid,
                self.s.tid,
                self.s.strand,
                f'Chr{self.s.chromosome}',
                self.s.biotype,
                f'tsl_{self.s.tsl}',
                self.s.counts,)
        else:
            self.header = self.header

    def level_changer(self):
        if self.id.startswith('sp'):
            if self.header.split('|')[9] == 'protein_coding':
                self.level = "L1"
            elif self.level == "L4":
                self.level = "L3"
            elif self.level == "L5":
                self.level = "L3"


    def str_to_seqrec(self):
        self.rec = SeqRecord(Seq(self.p.prot),self.header,description=self.s.gene_symbol)
        # print(self.rec)
        return self.rec

class ORFs:
    def __init__(self,
                 sequence: Sequences,
                 peptide: Peptide,
                 CanDB,
                 IsoDB):
        self.orfs = peptide.dict2

        self.header = peptide.header.split('|')
        self.header_orf = peptide.header
        self.level = 'ORF'
        self.gene_symbol = self.header[1]
        self.gid = self.header[2]
        self.tid = self.header[3]
        self.strand = self.header[4]
        self.chromosome = self.header[5]
        self.biotype = self.header[6]
        self.tsl = self.header[7]
        self.counts = self.header[8]
        self.id = '-'
        self.CanDB = CanDB
        self.IsoDB = IsoDB


    @staticmethod
    def string_fix(lst):
        return lst[0]

    def dict_parse(self):


        converted = pd.DataFrame.from_dict(self.orfs, orient='index')
        converted['phase'] = converted['phase'].apply(ORFs.string_fix)
        converted['len'] = converted['len'].apply(ORFs.string_fix)
        converted['start'] = converted['start'].apply(ORFs.string_fix)
        converted['stop'] = converted['stop'].apply(ORFs.string_fix)
        for i in converted.index:
            if i == '':
                converted = converted.drop(i)
        max_value = converted['len'].max()
        converted = converted[converted.len != max_value]
        self.converted = converted

    def write_header_loop(self, out_location, prefix):

        for i in self.converted.index:
            phase = self.converted.loc[i]['phase']
            length = self.converted.loc[i]['len']
            start = self.converted.loc[i]['start']
            stop = self.converted.loc[i]['stop']
            pep = i
            header = "{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|{10}|{11}|{12}".format(
                self.level,
                self.gene_symbol,
                self.gid,
                self.tid,
                self.strand,
                f'Chr{self.chromosome}',
                self.biotype,
                self.tsl,
                self.counts,
                f'Len={length}',
                f'Start_Site_{start}',
                f'Stop_Site_{stop}',
                f'{phase}',)


            self.canonical_aa = ''
            self.old = self.level
            with open(self.CanDB) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    if record.seq == pep:
                        self.level = 'Canonical'
                        self.id = record.id


            self.canonical_aa = ''
            self.old = self.level
            with open(self.IsoDB) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    if record.seq == pep:
                        self.id = record.id
            if self.id != '-':
                header = self.id + "{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|{10}|{11}|{12}".format(
                    self.level,
                    self.gene_symbol,
                    self.gid,
                    self.tid,
                    self.strand,
                    f'Chr{self.chromosome}',
                    self.biotype,
                    self.tsl,
                    self.counts,
                    f'Len={length}',
                    f'Start_Site_{start}',
                    f'Stop_Site_{stop}',
                    f'{phase}', )
            else:
                header = self.header_orf

            if self.id.startswith('sp'):
                if self.header.split('|')[9] == 'protein_coding':
                    self.level = "L1"
                elif self.level == "L4":
                    self.level = "L3"
                elif self.level == "L5":
                    self.level = "L3"

            rec = SeqRecord(Seq(pep), f'{header}', description=self.gene_symbol)

            lrf.prot_to_fasta(rec, out_location, prefix, "_altORFs")














