# -*- coding: utf-8 -*-
import argparse
import insilico_translation as ist
import longread_functions as lrf
import post_run_stats as prs
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from Bio import SeqIO
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Server # for retrieval of Uniprot IDs
import constants
import os
from time import sleep
from joblib import Parallel, delayed
import multiprocessing
import threading
import concurrent.futures

# Main function
def main():
    start = time.perf_counter()
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', help='path to gtf file', required=True)
    parser.add_argument('-f', '--fasta', help='path to fasta file', required=True)
    parser.add_argument('-p', '--file_prefix', help='filename prefix (not path)', required=True)
    parser.add_argument('-o', '--outpath', help='Output path. Note JCASTLR autonames files', required=True)
    args = parser.parse_args()
    # set constants
    gtf = args.gtf
    fasta = args.fasta
    prefix = args.file_prefix
    out_location = args.outpath
    g = ist.Gtf(gtf,'results/DGE/counts_matrix.tsv')
    g.read_cutoff('results')
    duplicate_count = 0
    with open(fasta) as f:
        n = 0
        for record in SeqIO.parse(f, 'fasta'):
            n += 1
    print(f'{n} transcripts to process.')

    # loop through each record to get info and translate
    num_cores = ((multiprocessing.cpu_count() - 2)*2)
    print(num_cores)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
        with open(fasta) as f:
            n1 = 0
            for record in SeqIO.parse(f, 'fasta'):
                n1 += 1
                lrf.progress_bar(n1, n, 50)
                executor.submit(paralell_me(record, g, out_location, prefix), record)
    sleep(7)
    prs.post_run_counts(out_location,prefix)
    print(f'{time.perf_counter() - start} seconds')

def paralell_me(record,g,out_location, prefix):
    # n1 += 1
    # lrf.progress_bar(n1, n, 50)
    r = record
    s = ist.Sequences(g, r)
    s.get_counts()
    if s.counts <= g.min_count:
        return 1
    s.subset_gtf()
    s.get_meta()
    # a = ist.Canonical_test(s)
    # a.get_canonical_aa()
    # a.make_header()
    # print(a.level, a.rid, a.biotype)
    if s.level == "Canonical":
        # a = ist.Canonical_test(s)
        # a.get_canonical_aa()
        # a.make_header()
        s.annotated_trancript_trim()
        p = ist.Peptide(s)
        prot_seq = p.annotated_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # seq = p.str_to_seqrec(prot_seq)
        # Canonical.append(seq)
    elif s.level == "L1":
        # a = ist.Canonical_test(s)
        # a.get_canonical_aa()
        # a.make_header()
        s.annotated_trancript_trim()
        p = ist.Peptide(s)
        prot_seq = p.annotated_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # print(prot_seq)
        # seq = p.str_to_seqrec(prot_seq)
        # L1.append(seq)
    elif s.level == "L2":
        # a = ist.Canonical_test(s)
        # a.get_canonical_aa()
        # a.make_header()
        s.annotated_trancript_trim()
        p = ist.Peptide(s)
        p1 = p.annotated_translate()
        p2 = p.multi_phase_translate()
        if len(p1) > len(p2):
            p.prot = p1
        elif len(p1) < len(p2):
            p.prot = p2
        else:
            p.prot = p1
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # seq = p.str_to_seqrec(prot_seq)
        # L2.append(seq)
    elif s.level == "L3":
        # a = ist.Canonical_test(s)
        # a.get_canonical_aa()
        # a.make_header()
        p = ist.Peptide(s)
        prot_seq = p.multi_phase_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # seq = p.str_to_seqrec(prot_seq)
        # L3.append(seq)
    elif s.level == "L4":
        # a = ist.Canonical_test(s)
        # a.get_canonical_aa()
        # a.make_header()
        p = ist.Peptide(s)
        prot_seq = p.multi_phase_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # seq = p.str_to_seqrec(prot_seq)
        # L4.append(seq)
    elif s.level == "L5":
        # a = ist.Canonical_test(s)
        # a.make_header()
        p = ist.Peptide(s)
        prot_seq = p.multi_phase_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        # seq = p.str_to_seqrec(prot_seq)
        # L5.append(seq)
    else:
        print("Orphan Read")

    ph = ist.Post_hoc_reassignment(s, p)
    ph.get_canonical_aa_uniprot_local('resources/DB/reviewed_canonical.fasta')
    ph.make_header()

    # post hoc change of level
    if ph.level == 'Canonical':
        seq = ph.str_to_seqrec()
        lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")

        # if ph.old == 'L1':
        #     # val = L1[-1]
        #     seq = ph.str_to_seqrec()
        #     # print("here   ",val)
        #     lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")
        # elif ph.old == 'L2':
        #     # val = L2[-1]
        #     seq = ph.str_to_seqrec()
        #     lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")
        # elif ph.old == 'L3':
        #     # val = L3[-1]
        #     seq = ph.str_to_seqrec()
        #     lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")
        # elif ph.old == 'L4':
        #     # val = L4[-1]
        #     seq = ph.str_to_seqrec()
        #     lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")
        # elif ph.old == 'L5':
        #     # val = L5[-1]
        #     seq = ph.str_to_seqrec()
        #     lrf.prot_to_fasta(seq, out_location, prefix, "_Canonical")
        # else:
        #     print("error post hoc")
    elif ph.level != 'Canonical':
        ph.get_aa_uniprot_local('resources/DB/reviewed_included_isoforms.fasta')
        ph.make_header()
        seq = ph.str_to_seqrec()
        if ph.level == 'L1':
            # val = L1[-1]
            # seq = ph.str_to_seqrec()
            # print("here   ",val)
            lrf.prot_to_fasta(seq, out_location, prefix, "_Level1")
        elif ph.level == 'L2':
            # val = L2[-1]
            # seq = ph.str_to_seqrec()
            lrf.prot_to_fasta(seq, out_location, prefix, "_Level2")
        elif ph.level == 'L3':
            # val = L3[-1]
            # seq = ph.str_to_seqrec()
            lrf.prot_to_fasta(seq, out_location, prefix, "_Level3")
        elif ph.level == 'L4':
            # val = L4[-1]
            # seq = ph.str_to_seqrec()
            lrf.prot_to_fasta(seq, out_location, prefix, "_Level4")
        elif ph.level == 'L5':
            # val = L5[-1]
            # seq = ph.str_to_seqrec()
            lrf.prot_to_fasta(seq, out_location, prefix, "_Level5")
        else:
            print("error post hoc")
    else:
        print("Error")


