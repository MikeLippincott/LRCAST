import Bio
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import argparse
import os

# # Parse Arguments
# parser = argparse.ArgumentParser()
# parser.add_argument('-c', '--can', help='.fasta file path', required=True)
# parser.add_argument('-j', '--jcast', help='Database name, .fasta', required=True)
# parser.add_argument('-n', '--name1', help='.fasta file path', required=True)
# parser.add_argument('-m', '--name2', help='Database name, .fasta', required=True)
# parser.add_argument('-p', '--prefix', help='Output path. Note JCASTLR autonames files', required=True)
# parser.add_argument('-s', '--suffix', help='Output path. Note JCASTLR autonames files', required=True)
# parser.add_argument('-o', '--outpath', help='Output path. Note JCASTLR autonames files', required=True)
# parser.add_argument('-r', '--read_cutoff', help='cut off of read count', required=True)
# args = parser.parse_args()


def fasta_merge(can,
                  jcast,
                  name1,
                  name2,
                  prefix,
                  suffix,
                output,
                cutoff):
    """
    :param output: output path
    :param suffix: string to be added to filename before extension
    :return: fasta file
    """


    filtered_lst = []
    with open(jcast) as jc:
        for record_jc in SeqIO.parse(jc, 'fasta'):
            if record_jc.id.startswith('sp'):
                # print(record_jc.id.split('|')[11])
                if int(record_jc.id.split('|')[11]) > int(cutoff):
                    filtered_lst.append(record_jc)
            elif int(record_jc.id.split('|')[8]) > int(cutoff):
                filtered_lst.append(record_jc)
                # print(record_jc.id.split('|')[8])


    outfile = os.path.join(output,  prefix + str(cutoff) + '_' + name1 + '_' + name2 + '_' + suffix + '.fasta')
    # if the file already exists, open it and amend that record.
    with open(outfile, 'a') as output_handle:
        for i in filtered_lst:
            SeqIO.write(i, output_handle, 'fasta')


    existing_records = {}
    for existing_record in SeqIO.parse(outfile, 'fasta'):
        existing_records[existing_record.seq] = None

    with open(can) as c:
        for record_can in SeqIO.parse(c, 'fasta'):
            # for jc_record in existing_records:
            #     if jc_record.seq == record_can.seq:
            if record_can.seq in existing_records:
                    # print('duplicate')
                pass
            else:
                with open(outfile, 'a') as output_handle:
                    SeqIO.write(record_can, output_handle, 'fasta')


def main():
    # for i in range(0,550,50):
    fasta_merge('resources/DB/reviewed_canonical.fasta',
                'aORFsJCASTLR_C_1_2_3_4_5.fasta',
                'canonical' ,
                'jcast_all' ,
                'aORFs' ,
                'all_levels',
                '.',
                0)
    fasta_merge('resources/DB/reviewed_included_isoforms.fasta',
                'aORFsJCASTLR_C_1_2_3_4_5.fasta',
                'canonical+isoform',
                'jcast_all',
                'aORFs',
                'all_levels',
                '.',
                0)

if __name__ == '__main__':
    main()




def DB_generate():

    can = 'resources/DB/reviewed_canonical.fasta'
    alt = 'resources/DB/reviewed_included_isoforms.fasta'

    can_count = 0
    alt_count = 0
    new_alt_count = 0

    # turn each fasta into a series
    can_lst = []
    with open(can) as c:
        for can_rec in SeqIO.parse(c, 'fasta'):
            can_count += 1
            can_lst.append(can_rec.id)
    can_series = pd.DataFrame({'x':can_lst})

    alt_lst = []
    with open(alt) as a:
        for alt_rec in SeqIO.parse(a, 'fasta'):
            alt_count += 1
            alt_lst.append(alt_rec.id)
    alt_series = pd.DataFrame({'x':alt_lst})

    # Merge on alt isoforms
    df = pd.merge(can_series,alt_series, how='outer',indicator=True)
    alt_ids = df[df['_merge'] == 'right_only']
    new_alt_count = len(alt_ids)


    print("looping records")

    write_lst = []
    with open(alt) as a:
        for alt_rec in SeqIO.parse(a, 'fasta'):
            for i in alt_ids['x']:
                if alt_rec.id == i:
                    write_lst.append(alt_rec)

    print('Writing to File')
    outfile = 'resources/DB/reviewed_alternative_isoforms.fasta'
    for i in write_lst:
        with open(outfile, 'a') as output_handle:
            seqrecord = i
            SeqIO.write(seqrecord, output_handle, 'fasta')

    print(f'{can_count} Canonical records')

    print(f'{alt_count} Canonical and Alternative Isoforms')

    print(f'{new_alt_count} Alternative Isoforms')



# i = 1
# while i < 65:
#     if i == 1:
#         i += 1
#         print(i)
#     else:
#         i = i*2
#         print(i)
#
# print(i)


