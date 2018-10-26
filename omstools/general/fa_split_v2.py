from __future__ import division
import click
from Bio import SeqIO
import pandas as pd
import envoy
import os
from omstools.utils.config import CLICK_CONTEXT_SETTINGS
from omstools.utils.config import MutuallyExclusiveOption
import math


def build_fa_index(fa):
    fa_index = '{fa}.fai'.format(fa=fa)
    if not os.path.isfile(fa_index):
        print 'Building fa index'
        envoy.run('samtools faidx {fa}'.format(fa=fa))
    fai_df = pd.read_table(fa_index, header=None)
    return len(fai_df)


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.argument(
    'input_fa',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.argument(
    'out_dir',
    type=click.Path(file_okay=False),
    required=True,
)
@click.option(
    '--genome',
    is_flag=True,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['file_number', 'fa_number'],
)
@click.option(
    '--fa_number',
    type=click.INT,
    help='fasta number in one file.',
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['file_number', 'genome'],
)
@click.option(
    '--file_number',
    type=click.INT,
    help='file number of split fastas.',
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['fa_number', 'genome']
)
def main(input_fa, out_dir, genome, fa_number,
         file_number):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    total_number = build_fa_index(input_fa)
    seq_per_file = 0
    file_count = 0
    fa_dict = dict()
    if fa_number:
        file_number = math.ceil(total_number / fa_number)
    if genome:
        file_number = total_number
    for n, seq_rd in enumerate(SeqIO.parse(input_fa, 'fasta')):
        m = int(n % file_number)
        if genome:
            file_name = os.path.join(out_dir,
                                     '{s}.fa'.format(s=seq_rd.id))
        else:
            file_name = os.path.join(out_dir,
                                     'split_{}.fa'.format(m))
        fa_dict.setdefault(file_name, []).append(seq_rd)
    for each_file in fa_dict:
        SeqIO.write(fa_dict[each_file], each_file, 'fasta')


if __name__ == '__main__':
    main()
