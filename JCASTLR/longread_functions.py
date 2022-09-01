# -*- coding: utf-8 -*-
import constants
import os
from Bio import SeqIO

# Functions

# Translate function nt -> aa
def translate(nt,
              phase):
    """
    :param nt: string of nucleotides
    :param phase: int [0, 1, 2] translating frame to use to search for ORFs
    :return: amino acid sequence

    function translates nucleotide sequence into peptide with no start codon info.
    Function works by finding ORF (Start codon), and translating from start site
    """
    # TODO: Find alternative ORFs and translate
    # Codon Lookup table
    code = constants.genetic_code
    pep = ''
    for i in range(phase, len(nt) - 2, 3):
        aa = code[nt[i:i + 3]]
        if len(pep) == 0:
            if aa == 'M':
                pep += aa
        elif len(pep) > 0:
            if aa == 'X':
                return pep
            else:
                pep += aa
    return pep

# uses tranlate function to tranlate nt in 3 different frames
def multi_phase_translate(nt):
    """
    :param nt: string of nucleotides
    :return: longest peptide
    """
    d = {}
    for j in [0,1,2]:
        a = translate(nt,j)
        d[len(a)] = a
    return d[max(d)]

# export function peptide sequence to fasta
def prot_to_fasta(seqrecord,
                  output,
                  prefix,
                  suffix):
    """
    :param seqrecord: BioSeq seqrecord record object
    :param output: output path
    :param suffix: string to be added to filename before extension
    :return: fasta file
    """
    outfile = os.path.join(output, prefix + 'JCASTLR' + suffix + '.fasta')
    # if the file already exists, open it and amend that record.
    existing_records = []
    if os.path.exists(outfile):
        for existing_record in SeqIO.parse(outfile, 'fasta'):
            existing_records.append(existing_record)
    else:
        with open(outfile, 'w') as f:
            SeqIO.write(seqrecord, f, 'fasta')

    for read_record in existing_records:
        if read_record.seq == seqrecord.seq:
            return True

    with open(outfile, 'a') as output_handle:
        SeqIO.write(seqrecord, output_handle, 'fasta')

# progress bar for run time estimation
def progress_bar(current,
                 total,
                 bar_length=20):
    """
    :param current: current n of iteration
    :param total: total len(list) of iteration
    :param bar_length: how long the bar shall appear in terminal
    :return: progress arrow
    """
    fraction = current / total
    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '
    ending = '\n' if current == total else '\r'
    print(f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)
