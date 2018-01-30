
import click
from crimson import fastqc
import re
import os
import pandas as pd
import glob
import sys


def get_sample_name(filename, names):
    match_length = 0
    name = ''
    for each in names:
        if re.match(each, filename):
            if len(each) > match_length:
                match_length = len(each)
                name = each
    return name


@click.command()
@click.option(
    '-d',
    '--qc_dir',
    type=click.Path(file_okay=False),
    required=True,
    help='fastqc directory.')
@click.option(
    '-n',
    '--names',
    type=click.Path(dir_okay=False),
    required=True,
    help='sample name file.')
@click.option(
    '-o',
    '--output',
    default='fastqc.summary.txt',
    type=click.Path(dir_okay=False),
    help='output summary file path.')
def main(qc_dir, names, output):
    qc_dict = dict()
    sample_names = [each.strip() for each in open(names) if each.strip()]
    fastqc_files = [each for each in os.listdir(
        qc_dir) if each.endswith('.zip')]
    if not fastqc_files:
        fastqc_files = glob.glob('{d}/*fastqc/fastqc_data.txt'.format(
            d=qc_dir))
        if not fastqc_files:
            sys.exit('Can not find fastqc result in [{d}].'.format(
                d=qc_dir))
    for each_file in fastqc_files:
        if each_file[-4:] == '.zip':
            each_file_path = os.path.join(qc_dir, each_file)
            sample_name = get_sample_name(each_file, sample_names)
        else:
            each_file_path = each_file
            each_file_dir = os.path.basename(os.path.dirname(each_file))
            sample_name = get_sample_name(each_file_dir, sample_names)
        each_qc = fastqc.parse(each_file_path)
        reads_num = each_qc['Basic Statistics']['contents']['Total Sequences']
        max_reads_len = max([int(str(each['Length']).split('-')[0]) for each in
                             each_qc['Sequence Length Distribution']['contents']])
        min_reads_len = min([int(str(each['Length']).split('-')[0]) for each in
                             each_qc['Sequence Length Distribution']['contents']])
        if min_reads_len == max_reads_len:
            reads_len = min_reads_len
        else:
            reads_len = '{mi}-{mx}'.format(mi=min_reads_len, mx=max_reads_len)
        gc_reads = each_qc['Basic Statistics']['contents']['%GC'] * reads_num
        base_num = sum([each['Count'] * int(str(each['Length']).split('-')[0])
                        for each in
                        each_qc['Sequence Length Distribution']['contents']])
        q20_reads = sum([each['Count'] for each in
                         each_qc['Per sequence quality scores']['contents']
                         if each['Quality'] >= 20])
        q30_reads = sum([each['Count'] for each in
                         each_qc['Per sequence quality scores']['contents']
                         if each['Quality'] >= 30])
        dup_per = 100 - \
            each_qc['Sequence Duplication Levels']['Total Deduplicated Percentage']
        dup_reads = dup_per * reads_num
        qc_dict.setdefault('Sample', []).append(sample_name)
        qc_dict.setdefault('Reads_number', []).append(reads_num)
        qc_dict.setdefault('read_len', []).append(reads_len)
        qc_dict.setdefault('Data_size', []).append(base_num)
        qc_dict.setdefault('gc', []).append(gc_reads)
        qc_dict.setdefault('q20_reads', []).append(q20_reads)
        qc_dict.setdefault('q30_reads', []).append(q30_reads)
        qc_dict.setdefault('dup_reads', []).append(dup_reads)
    qc_df = pd.DataFrame(qc_dict)
    summary_qc = qc_df.groupby(['Sample'])[
        'Reads_number', 'Data_size',
        'gc', 'q20_reads', 'q30_reads', 'dup_reads']
    sqc_df = summary_qc.sum()
    sqc_df.loc[:, 'GC(%)'] = sqc_df.loc[:, 'gc'] / \
        sqc_df.loc[:, 'Reads_number']
    sqc_df.loc[:, 'Q20(%)'] = sqc_df.loc[:, 'q20_reads'] / \
        sqc_df.loc[:, 'Reads_number'] * 100
    sqc_df.loc[:, 'Q30(%)'] = sqc_df.loc[:, 'q30_reads'] / \
        sqc_df.loc[:, 'Reads_number'] * 100
    sqc_df.loc[:, 'Duplication(%)'] = sqc_df.loc[:, 'dup_reads'] / \
        sqc_df.loc[:, 'Reads_number']
    sqc_df.loc[:, 'Reads_number(M)'] = sqc_df.loc[:, 'Reads_number'] / \
        (10 ** 6)
    sqc_df.loc[:, 'Data_size(G)'] = sqc_df.loc[:, 'Data_size'] / \
        (10 ** 9)
    sqc_df.loc[:, 'Reads_length(bp)'] = qc_df.groupby(['Sample'])[
        'read_len'].max()
    out_col = ['Reads_number(M)', 'Reads_length(bp)', 'Data_size(G)',
               'GC(%)', 'Q20(%)', 'Q30(%)', 'Duplication(%)']
    out_qc_df = sqc_df.loc[:, out_col]
    out_qc_df.to_csv(output, sep='\t', float_format='%.2f')


if __name__ == '__main__':
    main()
