import click
import pandas as pd
import os


def get_rna_matrix_inf(m_file):
    with open(m_file) as file_inf:
        flag = 0
        for eachline in file_inf:
            if 'PF_BASES' in eachline:
                flag = 1
                continue
            if flag:
                eachline_inf = eachline.strip().split()
                return eachline_inf


OUT_HEADER = ['Aligned_bases',
              '% mRNA bases',
              '% intronic bases',
              '% intergenic bases',
              '% Strand specific reads']


@click.command()
@click.option('-d', '--analysis_dir', type=click.Path(exists=True),
              help='rnaseq matrix output directory.', required=True)
@click.option('-s', '--sample_inf', type=click.Path(exists=True),
              help='sample information, first column is sample id.',
              required=True)
@click.option('-m', '--mapping_method', type=click.Choice(['star', 'hisat2']),
              help='mapping software.', default='star')
def main(analysis_dir, sample_inf, mapping_method):
    sample_df = pd.read_table(sample_inf, header=None, index_col=1)
    summary_file = os.path.join(analysis_dir, 'rna_matrix.summary.txt')
    summary_file_inf = open(summary_file, 'w')
    header = ['{e} [{m}]'.format(e=each, m=mapping_method) for
              each in OUT_HEADER]
    summary_file_inf.write('Sample\t{h}\n'.format(h='\t'.join(header)))
    for each_sample in sample_df.index:
        each_matrix = os.path.join(
            analysis_dir, '{s}.RNA_Metrics'.format(s=each_sample))
        each_m_inf = get_rna_matrix_inf(each_matrix)
        summary_file_inf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
            each_sample, each_m_inf[1], each_m_inf[20], each_m_inf[18],
            each_m_inf[19], each_m_inf[22]
        ))
    summary_file_inf.close()


if __name__ == '__main__':
    main()
