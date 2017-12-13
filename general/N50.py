from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '--samtools_index',
    help='samtools index of genome',
    required=True)
parser.add_argument(
    '--cutoff',
    help='shortest scaffold/contig length.',
    default=1000,
    type=int)
args = parser.parse_args()

samtools_index = args.samtools_index

length_list = []

with open(samtools_index) as s_index_inf:
    for eachline in s_index_inf:
        eachline_inf = eachline.strip().split()
        chr_len = int(eachline_inf[1])
        if chr_len > args.cutoff:
            length_list.append(chr_len)


length_list.sort()
total_len = sum(length_list)
add_len = 0
for each_len in length_list:
    add_len += each_len
    if add_len > 0.5 * total_len:
        print('N50: %d' % each_len)
        break
