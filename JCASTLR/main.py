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
import multiprocessing as mp
from tqdm import tqdm

# Main function
def main():
    start = time.perf_counter()
    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", '--model', action='store_true')
    parser.add_argument('-r', '--read_cutoff',help='read cutoff value', required=False)
    parser.add_argument('-a', '--altORFs', action='store_true')
    parser.add_argument('-g', '--gtf', help='path to gtf file', required=True)
    parser.add_argument('-f', '--fasta', help='path to fasta file', required=True)
    parser.add_argument('-p', '--file_prefix', help='filename prefix (not path)', required=True)
    parser.add_argument('-o', '--outpath', help='Output path. Note JCASTLR autonames files', required=True)
    args = parser.parse_args()
    # set constants

    gtf = args.gtf
    global altORFs
    altORFs = args.altORFs
    # global it
    # it = 0
    fasta = args.fasta
    prefix = args.file_prefix
    out_location = args.outpath

    with open(fasta) as f:
        n = 0
        for record in SeqIO.parse(f, 'fasta'):
            n += 1
        print(f'{n} transcripts to process.')
    g = ist.Gtf(gtf,'results/DGE/counts_matrix.tsv')
    if args.model:
        g.read_cutoff('results')
        print("model")
    else:
        g.min_count = int(args.read_cutoff)
    print(g.min_count)

    duplicate_count = 0

    pool = mp.Pool(mp.cpu_count())
    n1 = 0
    with open(fasta) as f:
        # for record in SeqIO.parse(f, 'fasta'):
            # n1 += 1
            # lrf.progress_bar(n1, n, 50)
        results = pool.starmap_async(paralell_me,
                                     [(record,
                                            g,
                                            out_location,
                                            prefix) for record in SeqIO.parse(f, 'fasta')]).get()
        tqdm(total=n)

    pool.close()
    pool.join()

    # with open(fasta) as f:
    #     n1 = 0
    #     for record in SeqIO.parse(f, 'fasta'):
    #         n1 += 1
    #         lrf.progress_bar(n1, n, 50)
              # paralell_me(record, g, out_location, prefix)

    sleep(7)
    print('Starting Post Run Analysis')
    prs.post_run_counts(out_location,prefix,altORFs)
    print(f'{time.perf_counter() - start} seconds')


def paralell_me(record,g,out_location, prefix):
    # n1 += 1
    # lrf.progress_bar(n1, n, 50)
    r = record
    s = ist.Sequences(g, r)
    s.get_counts()
    # print(s.counts)
    # if int(s.counts) <= int(g.min_count):
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
        p1 = p.annotated_translate()
        p2 = p.multi_phase_translate()
        if len(p1) > len(p2):
            p.prot = p1
        elif len(p1) < len(p2):
            p.prot = p2
        else:
            p.prot = p1

        orfs = ist.ORFs(s,p,'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
        # seq = p.str_to_seqrec(prot_seq)
        # Canonical.append(seq)
    elif s.level == "L1":
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
        orfs = ist.ORFs(s, p, 'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
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
        orfs = ist.ORFs(s, p, 'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
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
        orfs = ist.ORFs(s, p, 'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
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
        orfs = ist.ORFs(s, p, 'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
        # seq = p.str_to_seqrec(prot_seq)
        # L4.append(seq)
    elif s.level == "L5":
        # a = ist.Canonical_test(s)
        # a.make_header()
        p = ist.Peptide(s)
        prot_seq = p.multi_phase_translate()
        # p.get_canonical_aa_uniprot_local()
        p.make_header()
        orfs = ist.ORFs(s, p, 'resources/DB/reviewed_canonical.fasta',
                               'resources/DB/reviewed_alternative_isoforms.fasta')
        orfs.dict_parse()
        orfs.write_header_loop(out_location, prefix)
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

    elif ph.level != 'Canonical':
        ph.get_aa_uniprot_local('resources/DB/reviewed_alternative_isoforms.fasta')
        ph.make_header()
        ph.level_changer()
        # if ph.old != ph.level:
            # print(ph.level,ph.s.level)

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
    # it += 1
    # print(it)
    print("completed")

