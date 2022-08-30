# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Server # for retrieval of Uniprot IDs
import constants
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

"""
Run this module with input files from the output of the FLAIR pipeline
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
        print("LOADING GTF FILES INTO PANDAS DATAFRAMES.")
        # FLAIR output gtf
        self.gtf = read_gtf(gtf_loc)
        # Ensembl gtf
        self.gtf_file = read_gtf("resources/genome/Homo_sapiens.GRCh38.107.gtf")

# Sequence object
class Sequences(object):
    def __init__(self,
                 annotation_df: Gtf,
                 record):
        self.gtf = annotation_df.gtf
        self.gtf_file = annotation_df.gtf_file
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
        self.frame = self.gtf0['frame'].unique()[0]
        if len(self.frame) == 1:
            self.frame = int(self.frame[0])
        else:
            self.frame == None
        # return self.id, self.strand, self.frame, self.transcript

    # retrieve transcript or gene meta data
    def get_meta(self):

        """
        Function defines each record's level then retrieves transcript or gene meta data depending on which ensembl id
        is available (ensembl transcript id has priority over ensembl gene id)
        Metadata to retieve:
        level
        biotyep
        gene symbol
        gene name
        chromosome
        transcript support level (tsl)
        uniprot id
        :return: No return
        """
        print(self.rid)
        if 'ENST' in self.rid:
            enst = gget.info(self.tid)
            gtf0 = self.gtf_file.query(f'transcript_id == "{self.tid}"').query('transcript_biotype == "protein_coding" ')
            if len(gtf0) > 0:
                self.level = "L1"
                self.biotype = "protein_coding"
            elif len(gtf0) == 0:
                self.biotype = np.unique(self.gtf_file.query(f'transcript_id == "{self.tid}"')['transcript_biotype'])
                self.level = "L3"
            else:
                print("Check Meta Data method in Sequences Class")

            sleep(0.1)

            try:
                if pd.isnull(enst):
                    self.gene_name = '-'
            except ValueError:
                try:
                    if pd.isnull(enst['protein_names'].to_list()):
                        self.gene_name = "-"
                    else:
                        self.gene_name = enst['protein_names'].to_list()[0]
                except ValueError:
                    # gene_name = "-"
                    self.gene_name = enst['protein_names'].to_list()[0][0]

            if self.gene_name == '-':
                enst = gget.info(self.gid)
                try:
                    if pd.isnull(enst):
                        self.gene_name = '-'
                except ValueError:
                    try:
                        if pd.isnull(enst['protein_names'].to_list()):
                            self.gene_name = "-"
                        else:
                            self.gene_name = enst['protein_names'].to_list()[0]
                    except ValueError:
                        # gene_name = "-"
                        self.gene_name = enst['protein_names'].to_list()[0][0]

            print(self.gene_name)
            self.gene_symbol = enst['ensembl_gene_name'].to_list()[0]
            self.chromosome = self.gtf0['seqname'].to_list()[0]
            if len(gtf0['transcript_support_level']) == 0:
                self.tsl = "-"
            else:
                self.tsl = gtf0['transcript_support_level'].to_list()[0]
            self.uniprot = enst['uniprot_id'].to_list()[0]
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        elif 'ENSG' in self.rid:
            gtf0 = self.gtf_file.query(f'gene_id == "{self.gid}"').query('transcript_biotype == "protein_coding" ')
            if len(gtf0) > 0:
                self.level = "L2"
                self.biotype = "protein_coding"
            elif len(gtf0) == 0:
                self.biotype = np.unique(self.gtf_file.query(f'gene_id == "{self.gid}"')['transcript_biotype'])
                self.level = "L4"
            else:
                print("Check Meta Data method in Sequences Class")

            # print(self.gid)
            ensg = gget.info(self.gid)
            sleep(0.1)
            # self.gene_name = ensg['protein_names'].to_list()
            # print(self.gene_name[0])
            # print(ensg['protein_names'])
            try:
                if pd.isnull(ensg):
                    self.gene_name = '-'
            except ValueError:
                try:
                    if pd.isnull(ensg['protein_names'].to_list()):
                        self.gene_name = "-"
                    else:
                        self.gene_name = ensg['protein_names'].to_list()[0]
                except ValueError:
                    # gene_name = "-"
                    self.gene_name = ensg['protein_names'].to_list()[0][0]

            print(self.gene_name)
            self.gene_symbol = ensg['ensembl_gene_name'].to_list()[0]
            self.chromosome = self.gtf0['seqname'].to_list()[0]
            self.uniprot = ensg['uniprot_id'].to_list()[0]
            if len(gtf0['transcript_support_level']) == 0 :
                self.tsl = "-"
            else:
                self.tsl = gtf0['transcript_support_level'].to_list()[0]
            # return self.level, self.gene_name, self.gene_symbol, self.uniprot, self.chromosome

        elif not "ENSG" in self.rid:
            self.level = "L5"
            self.biotype = "L5"
            self.gene_name = "L5"
            self.gene_symbol = "L5"
            self.chromosome = self.gtf0['seqname'].to_list()[0]
            self.uniprot = "L5"
            self.tsl = "L5"

        else:
            print('LRCAST: RNA Fasta Header Read Error!')

        if pd.isnull(self.gene_name):
            self.gene_name = '-'

    def annotated_trancript_trim(self):
        """
        :return: trimmed sequence for insilico translation
        """
        gtf0 = self.gtf_file.query(f'transcript_id == "{self.tid}"').query('transcript_biotype == "protein_coding" ')
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
        if len(TSS) > 0:
            y1 = TSS['start'].to_list()[0]
            y2 = TSS['end'].to_list()[0]
        elif len(TSS) == 0:
            print("No Annotated Start Codon")
            return None
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
        print(self.a_seq)
        # print(start, end)
        # print(e)
        # print(SS)
        # print(t_len)
        # print(adjustedSS)
        return self.a_seq

uniprot_max_retries = 10 # Defined for number of retries for uinprot address in Canonical_tests
class Canonical_test:
    def __init__(self,
                 sequence: Sequences):
        self.gtf_canonical_transcript = None
        self.gtf_alternative_transcripts = None
        self.s = sequence
        self.logger = logging.getLogger('jcast.seq')
        # self.canonical_aa = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

    def get_canonical_aa(self):
        """
        :return: True
        """

        self.canonical_aa = self.get_canonical_aa_uniprot(reviewed='true')
        if len(self.canonical_aa[:]) == 0:
            self.canonical_aa = self.get_canonical_aa_uniprot(reviewed='false')

        # print(str(self.canonical_aa))

        # if str(self.canonical_aa.seq) == '':
        #     self.canonical = 'Iso'
        # else:
        #     self.canonical = 'Canonical'
        # self.s.header += f'{canonical}|'
        # print(self.canonical_aa.id,self.canonical_aa.seq,self.canonical_aa.description)

        return True

    def get_canonical_aa_uniprot(self,
                                 reviewed='true',
                                 ) -> SeqRecord:

        """
        get the canonical sequences from Uniprot
        :param reviewed: get only reviewed sequence
        :return: canonical aa seqrecord
        """

        record = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

        # cache retrieved sequences if available, otherwise retrieve from Ensembl
        cache = self.s.gid

        # Create cache folder if not exists
        if not os.path.exists('cache'):
            os.makedirs('cache')

        con = sq.connect(os.path.join('cache', 'uniprot-cache.db'))
        cur = con.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')

        cur.execute('''SELECT id, seq FROM sequences WHERE id=:cache''',
                    {'cache': cache})
        read_fasta = cur.fetchone()

        if read_fasta:
            record = list(SeqIO.parse(StringIO(read_fasta[1]), 'fasta'))[0]
            self.logger.info('Locally cached sequence retrieved for {0}.'.format(cache))
            con.close()

        else:

            server = 'https://www.ebi.ac.uk'

            # Uniprot does not support transcript version
            ext = '/proteins/api/proteins/Ensembl:' + re.sub('\\..*', '', self.s.gid) + \
                  '?offset=0&size=1&reviewed=' + reviewed + '&isoform=0'

            self.logger.info('Sequence not cached locally. Attempting to get from Uniprot: {0}'.format(
                server + ext
            ))

            retries = Retry(total=uniprot_max_retries,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            rqs = rq.Session()
            rqs.mount('https://', HTTPAdapter(max_retries=retries))
            ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

            if not ret.ok:
                self.logger.warning('Failed retrieval for {0} after {1} retries.'.format(self.s.gid,
                                                                                         uniprot_max_retries)
                                    )
                con.close()  # TODO: change to with statement

            if ret.status_code == 200 and ret.text != '':
                record = list(SeqIO.parse(StringIO(ret.text), 'fasta'))[0]

                # TODO: Catch sqlite3.OperationalError if writing fails, such as in a remote volume.
                cur.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''',
                            {'id': cache, 'seq': record.format('fasta')})
                con.commit()
                con.close()
                self.logger.info('Sequence retrieved from Uniprot and written into local cache.')

            elif ret.status_code == 200 and ret.text == '':
                self.logger.info('Retrieved empty fasta from Ensembl for {0}'.format(self.s.gid))
                # TODO: A known issue where Uniprot sequences do not exist for some Ensembl genes
                # TODO: Rather than fix this we will simply use the GTF in future versions.
                con.close()

            elif ret.status_code != 200:
                self.logger.warning('Retrieval of protein sequence failed.')
                con.close()

            else:
                con.close()
                pass

        return record

    def make_header(self):
        if str(self.canonical_aa.seq) == '':
            self.header = "{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|".format(
                f'LRCAST{self.s.level}',
                self.s.uniprot,
                self.s.gene_symbol,
                self.s.gid,
                self.s.tid,
                self.s.strand,
                f'Chr{self.s.chromosome}',
                self.s.biotype,
                self.s.tsl,
                self.s.level, )
        else:
            self.header = self.canonical_aa.id + "|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|".format(
                self.s.uniprot,
                self.s.gene_symbol,
                self.s.gid,
                self.s.tid,
                self.s.strand,
                f'Chr{self.s.chromosome}',
                self.s.biotype,
                self.s.tsl,
                self.s.level, )

        # return self.header

# Peptide Seq object class
class Peptide(object):
    def __init__(self,
                 sequence: Sequences,
                 canonical: Canonical_test):
        # sequence.subset_gtf()
        # sequence.get_meta()
        # sequence.annotated_trancript_trim()
        self.sequence = sequence
        canonical.make_header()
        self.seq = str(sequence.seq)
        self.strand = sequence.strand
        self.frame = sequence.frame
        self.header = canonical.header
        self.gene_name = sequence.gene_name

    def annotated_translate(self):
        """
        used the trimmed sequence to translate from start codon
        :return: peptide from annotation coordinates
        """
        code = constants.genetic_code
        pep = ''
        for i in range(len(self.seq) - 2, 3):
            aa = code[self.sequence.a_seq[i:i + 3]]
            if aa == 'X':
                return pep
            else:
                pep += aa
        return pep

    @staticmethod
    def start_find_translate(nt, phase):
        """
        :nt: input nucleotide sequence
        :phase: integer 0-2 which phase to translate in
        :return: longest peptide
        """
        code = constants.genetic_code
        pep = ''
        for i in range(phase, len(nt) - 2, 3):
            aa = code[nt[i:i + 3]]
            if len(pep) == 0:
                if aa == 'M':
                    pep += aa
                    # print(i)
            elif len(pep) > 0:
                if aa == 'X':
                    return pep
                else:
                    pep += aa
        # print(len(pep),", ",pep)
        return pep

    # do a multiphase translation using the translate function to find ORFs
    def multi_phase_translate(self):
        """
        :return: longest peptide
        """
        d = {}
        for j in [0, 1, 2]:
            a = Peptide.start_find_translate(self.seq, j)
            d[len(a)] = a
            self.prot = d[max(d)]
        return self.prot

    # TODO: write method to discern which translate method to use based upon value of seq from sequences class.
    # TODO: also include a report to compare methods of translation?
    def which_translate(self):
        print("not Complete")

    # returns longest translated aa to a BioSeq Record object
    def str_to_seqrec(self, prot_seq):
        print(self.gene_name)
        self.rec = SeqRecord(Seq(prot_seq),f'{self.header}',description=self.gene_name)
        return self.rec
