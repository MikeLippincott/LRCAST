import Bio
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

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

