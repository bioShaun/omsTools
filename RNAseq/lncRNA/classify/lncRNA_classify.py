import sys
from HTSeq import GFF_Reader
import os
import click
import pandas as pd
try:
    import envoy
except ImportError:
    import subprocess


CURRENT_DIR = os.getcwd()
VARS = dir()


def _oms_lncRNA_feelnc(mrna, lncrna, dis=5000):
    '''
    run FEElnc_classifier.pl to get lncRNA position relative to mRNA
    in order to classify them.
    '''
    cmd_line = 'FEELnc_classifier.pl -i {l} -a {m} -w {d} -m {d}'.format(
        l=lncrna, m=mrna, d=dis
    )
    if 'envoy' in VARS:
        r = envoy.run(cmd_line)
        output = r.std_out
    else:
        cmd_list = cmd_line.split()
        p = subprocess.Popen(cmd_list,
                             universal_newlines=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             bufsize=0,
                             cwd=CURRENT_DIR)
        output = p.communicate()[0]
    return output


def _oms_lncRNA_classify(feelnc_prd, lnc_gtf, method='Luo'):
    '''
    add lncRNA class to FEElnc_classifier.pl output
    '''
    feelnc_df = pd.read_table(feelnc_prd, index_col=2)
    feelnc_best_df = feelnc_df[feelnc_df.isBest == 1]
    lnc_class_df = feelnc_best_df.loc[:, feelnc_best_df.columns[1:]]
    lnc_class_df.loc[:, 'classification'] = None
    for eachline in GFF_Reader(lnc_gtf):
        tr_id = eachline.attr['transcript_id']
        gene_id = eachline.attr['gene_id']
        if tr_id in lnc_class_df.index:
            if not lnc_class_df.loc[tr_id, 'classification'] is None:
                continue
            if method == 'Luo':
                dirt, ltype, dis, subtype, loc = lnc_class_df.loc[tr_id][3:-1]
                lnc_class = _get_luo_code(dirt, subtype, loc)
                lnc_class_df.loc[tr_id, 'classification'] = lnc_class
            else:
                sys.exit('undefined classification method.')
        else:
            class_detail = ['--' for each in lnc_class_df.columns]
            class_detail[0] = gene_id
            class_detail[-1] = 'lincRNA'
            lnc_class_df.loc[tr_id] = class_detail
    return lnc_class_df


def _get_luo_code(*fee_loc):

    genic_dict = {
        'antisense': 'X',
        'sense': 'S',
        'overlapping': 'P',
        'nested': 'I',
        'containing': 'O',
        'divergent': 'H',
        'convergent': 'T',
        'same_strand': '',
        'downstream': 'D',
        'upstream': 'U',
        'intronic': 'I',
        'exonic': 'E',
    }

    code = [genic_dict[each] for each in fee_loc]
    code = ''.join(code)

    modify_dict = {
        'XHU': 'XH',
        'XTD': 'XT',
    }

    return modify_dict.get(code, code)


@click.command()
@click.option(
    '-m',
    '--mrna',
    type=click.Path(dir_okay=False, exists=True),
    help='mRNA gtf file.',
    required=True)
@click.option(
    '-l',
    '--lncrna',
    type=click.Path(dir_okay=False, exists=True),
    help='lncRNA gtf file.',
    required=True
)
@click.option(
    '-d',
    '--distance',
    type=click.INT,
    help='Size of the window around the lncRNA to compute classification.',
    default=5000
)
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    default=CURRENT_DIR,
    help='Directory to store classification output.'
)
def main(mrna, lncrna, out_dir, distance):
    '''
    classify lncRNA according to it's position relative to mRNA
    '''
    feelnc_out = os.path.join(out_dir, 'feelnc.classification.txt')
    if not os.path.exists(feelnc_out):
        output = _oms_lncRNA_feelnc(mrna, lncrna, distance)
        with open(feelnc_out, 'w') as out_inf:
            out_inf.write(output)
    lnc_class_df = _oms_lncRNA_classify(feelnc_out, lncrna)
    lnc_class_out = os.path.join(out_dir, 'lncRNA.classification.txt')
    lnc_class_df.to_csv(lnc_class_out, sep='\t')


if __name__ == '__main__':
    main()
