import os
import Bio
from Bio import SeqIO
import main




def post_run_counts(out_location, prefix, orf_var):
    fc = os.path.join(out_location, prefix + 'JCASTLR' + '_Canonical' + '.fasta')
    f1 = os.path.join(out_location, prefix + 'JCASTLR' + '_Level1' + '.fasta')
    f2 = os.path.join(out_location, prefix + 'JCASTLR' + '_Level2' + '.fasta')
    f3 = os.path.join(out_location, prefix + 'JCASTLR' + '_Level3' + '.fasta')
    f4 = os.path.join(out_location, prefix + 'JCASTLR' + '_Level4' + '.fasta')
    f5 = os.path.join(out_location, prefix + 'JCASTLR' + '_Level5' + '.fasta')
    dup = os.path.join(out_location, prefix + 'JCASTLR' + '_Duplicates' + '.fasta')
    aORF = os.path.join(out_location, prefix + 'JCASTLR' + '_altORFs' + '.fasta')



    with open(fc) as f:
        c = 0
        for record in SeqIO.parse(f, 'fasta'):
            c += 1
        print(f'{c} Canonical Isoforms')

    with open(f1) as f:
        l1 = 0
        for record in SeqIO.parse(f, 'fasta'):
            l1 += 1
        print(f'{l1} Level 1 Isoforms')

    with open(f2) as f:
        l2 = 0
        for record in SeqIO.parse(f, 'fasta'):
            l2 += 1
        print(f'{l2} Level 2 Isoforms')

    with open(f3) as f:
        l3 = 0
        for record in SeqIO.parse(f, 'fasta'):
            l3 += 1
        print(f'{l3} Level 3 Isoforms')

    with open(f4) as f:
        l4 = 0
        for record in SeqIO.parse(f, 'fasta'):
            l4 += 1
        print(f'{l4} Level 4 Isoforms')

    with open(f5) as f:
        l5 = 0
        for record in SeqIO.parse(f, 'fasta'):
            l5 += 1
        print(f'{l5} Level 5 Isoforms')

    with open(dup) as f:
        d = 0
        for record in SeqIO.parse(f, 'fasta'):
            d += 1
        print(f'{d} Duplicate Isoforms')

    if orf_var:

        with open(aORF) as f:
            a = 0
            for record in SeqIO.parse(f, 'fasta'):
                a += 1
            print(f'{a} altORF Isoforms')

        total = c + l1 + l2 + l3 + l4 + l5 + d + a
        print(f'{total} total Catorgized transcripts')

        unique_total = c + l1 + l2 + l3 + l4 + l5 + a
        print(f'{unique_total} total unique Catorgized transcripts')

    else:

        total = c + l1 + l2 + l3 + l4 + l5 + d
        print(f'{total} total Catorgized transcripts')

        unique_total = c + l1 + l2 + l3 + l4 + l5
        print(f'{unique_total} total unique Catorgized transcripts')