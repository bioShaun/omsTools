from __future__ import print_function
import pandas as pd
import click


@click.command()
@click.argument('chr_size', metavar='<chr_size>')
@click.argument('gtf', metavar='<gtf>')
def main(chr_size, gtf):
    '''
    taco output some transcripts whose ends exceeding the chromosome
    length, this script will find these transcripts and modify their
    end to chromosome end.
    '''
    chr_size_df = pd.read_table(chr_size, header=None, index_col=0)
    with open(gtf) as gtf_inf:
        for eachline in gtf_inf:
            eachline_inf = eachline.strip().split('\t')
            chrom = eachline_inf[0]
            end = int(eachline_inf[4])
            chrom_size = chr_size_df.loc[chrom][1]
            if end > chrom_size:
                eachline_inf[4] = str(chrom_size)
            out_line = '\t'.join(eachline_inf)
            print('{out_line}'.format(**locals()))


if __name__ == '__main__':
    main()
