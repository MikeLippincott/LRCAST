# -*- coding: utf-8 -*-
import argparse
import insilico_translation as ist
import longread_functions as lrf
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Server # for retrieval of Uniprot IDs
import constants
import os
from time import sleep

# Main function
def main():
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', help='path to gtf file', required=True)
    parser.add_argument('-b', '--bed', help='path to bed file', required=True)
    parser.add_argument('-f', '--fasta', help='path to fasta file', required=True)
    parser.add_argument('-p', '--file_prefix', help='filename prefix (not path)', required=True)
    parser.add_argument('-o', '--outpath', help='Output path. Note JCASTLR autonames files', required=True)
    args = parser.parse_args()
    # set constants
    gtf =  args.gtf
    bed = args.bed
    fasta = args.fasta
    prefix = args.file_prefix
    out_location = args.outpath
    g = ist.Gtf(gtf)
    b = ist.Bed(bed)
    m = ist.LoadMart()
    with open(fasta) as f:
        canon = []
        level1 = []
        level2 = []
        n = 0
        n1 = 0
        for line in f:
            if line.startswith(">"):
                n += 1
    print(fasta)
    with open(fasta) as f:
        for record in SeqIO.parse(f, 'fasta'):
            n1 += 1
            lrf.progress_bar(n1, n, 50)
            r = record
            a = ist.Sequences(g, m, r)
            a.subset_gtf()
            a.get_meta()
            a.make_header()
            if a.level == "canonical":
                p = ist.Peptide(a)
                p.multi_phase_translate()
                seq = p.str_to_seqrec()
                canon.append(seq)
            elif a.level == "Level 1":
                p = ist.Peptide(a)
                p.multi_phase_translate()
                seq = p.str_to_seqrec()
                level1.append(seq)
            elif a.level == "Level 2":
                p = ist.Peptide(a)
                p.multi_phase_translate()
                seq = p.str_to_seqrec()
                level2.append(seq)
            else:
                print("a.id")
        for i in canon:
            lrf.prot_to_fasta(i, out_location, prefix,"canonical")
        print(f'{len(canon)} canonical Isoforms')
        for i in level1:
            lrf.prot_to_fasta(i, out_location,prefix, "level1")
        print(f'{len(level1)} Level 1 Isoforms')
        for i in level2:
            lrf.prot_to_fasta(i, out_location, prefix,"level2")
        print(f'{len(level2)} Level 2 Isoforms')






