from __future__ import division
import click
import pandas as pd
import os
import glob


STAR_COLUMNS = ['Average input read length',
                'Number of input reads',
                'Uniquely mapped reads number',
                'Uniquely mapped reads %',
                'Number of reads mapped to multiple loci',
                '% of reads mapped to multiple loci']

# STAR_OUT_COL = ['Mean read length (bp)',
#                 'No. of reads (millions)',
#                 'Total no. of read alignments (millions)',
#                 '% of reads aligned',
#                 'No. of non-unique alignments (millions)',
#                 '% of non-unique alignments']


STAR_OUT_COL = ['Mean read length (bp)',
                'No. of reads (millions)',
                '% of reads aligned',
                '% of non-unique alignments']


def read_star_mapping_log(star_log_file):
    star_df = pd.read_table(star_log_file, header=None, sep='|', index_col=0)
    star_df = star_df.dropna()
    star_df = star_df.ix[4:]
    star_df.loc[:, 1] = [each.strip() for each in star_df.loc[:, 1]]
    return star_df


def add_percent(num):
    return '{n:.2f}%'.format(n=num)


@click.command()
@click.option('-s', '--sample_inf', type=click.Path(exists=True),
              help='sample information, second column is sample id.',
              required=True,)
@click.option('-d', '--mapping_dir', type=click.Path(exists=True),
              help='mapping directory.', required=True)
@click.option('-m', '--mapping_method', type=click.Choice(['star', 'hisat2']),
              help='mapping software.', default='star')
def main(sample_inf, mapping_dir, mapping_method):
    sample_df = pd.read_table(sample_inf, header=None, index_col=1)
    if mapping_method == 'star':
        out_file = os.path.join(mapping_dir, 'star_mapping.summary.txt')
        star_log_file_list = [os.path.join(
            mapping_dir, each_sample, 'Log.final.out')
            for each_sample in sample_df.index]
        star_log_df_list = map(read_star_mapping_log, star_log_file_list)
        star_log_df = pd.concat(star_log_df_list, axis=1)
        star_log_df.columns = sample_df.index
        star_log_out_df = star_log_df.T
        star_log_out_df.index.name = 'Sample'
        star_log_out_df.columns = [each.strip()
                                   for each in star_log_out_df.columns]
        out_data_tmp = star_log_out_df.loc[:, STAR_COLUMNS]
        out_data_tmp.loc[:, 'Number of input reads'] = map(
            int, star_log_out_df.loc[:, 'Number of input reads'])
        out_data_tmp.loc[:, 'Uniquely mapped reads number'] = map(
            int, out_data_tmp.loc[:, 'Uniquely mapped reads number'])
        out_data_tmp.loc[:, 'Number of reads mapped to multiple loci'] = map(
            int, out_data_tmp.loc[:, 'Number of reads mapped to multiple loci'])
        out_data_tmp.loc[:, 'Total no. of read alignments (millions)'] = (
            out_data_tmp.loc[:, 'Uniquely mapped reads number'] +
            out_data_tmp.loc[:, 'Number of reads mapped to multiple loci'])
        out_data_tmp.loc[:, 'Mean read length (bp)'] = out_data_tmp.loc[:, 'Average input read length']
        out_data_tmp.loc[:, 'No. of reads (millions)'] = out_data_tmp.loc[:, 'Number of input reads'] / 10**6
        out_data_tmp.loc[:, '% of reads aligned'] = 100 * out_data_tmp.loc[:, 'Total no. of read alignments (millions)'] / out_data_tmp.loc[:, 'Number of input reads']
        out_data_tmp.loc[:, 'Total no. of read alignments (millions)'] = out_data_tmp.loc[:, 'Total no. of read alignments (millions)'] / 10**6
        out_data_tmp.loc[:, '% of reads aligned'] = map(add_percent, out_data_tmp.loc[:, '% of reads aligned'])
        out_data_tmp.loc[:, 'No. of non-unique alignments (millions)'] = out_data_tmp.loc[:, 'Number of reads mapped to multiple loci'] / 10**6
        out_data_tmp.loc[:, '% of non-unique alignments'] = out_data_tmp.loc[:, '% of reads mapped to multiple loci']
        out_data_tmp.loc[:, 'No. of reads (millions)'] = out_data_tmp.loc[:, 'No. of reads (millions)'].map('{:,.2f}'.format)
        out_data_tmp.loc[:, 'Total no. of read alignments (millions)'] = out_data_tmp.loc[:, 'Total no. of read alignments (millions)'].map('{:,.2f}'.format)
        out_data_tmp.loc[:, 'No. of non-unique alignments (millions)'] = out_data_tmp.loc[:, 'No. of non-unique alignments (millions)'].map('{:,.2f}'.format)
        out_data_tmp.columns = ['{e} [star]'.format(e=each) for each in out_data_tmp.columns]
        star_out_col = ['{e} [star]'.format(e=each) for each in STAR_OUT_COL]
        out_data_tmp.to_csv(out_file, columns=star_out_col, sep='\t')
    if mapping_method == 'hisat2':
        script_dir = os.path.join(mapping_dir, 'script')
        out_file = os.path.join(mapping_dir, 'hisat2_mapping.summary.txt')
        out_file_inf = open(out_file, 'w')
        # out_file_inf.write('Sample\tTotal no. of read alignments (millions) [hisat2]\t% of reads aligned [hisat2]\tNo. of non-unique alignments (millions) [hisat2]\t% of non-unique alignments [hisat2]\n')
        out_file_inf.write('Sample\t% of reads aligned [hisat2]\t% of non-unique alignments [hisat2]\n')
        log_files = glob.glob('{d}/*log.*'.format(d=script_dir))
        for each_log in log_files:
            with open(each_log) as log_inf:
                for eachline in log_inf:
                    if 'were paired' in eachline:
                        total = 2 * int(eachline.strip().split()[0])
                    if 'aligned concordantly exactly 1 time' in eachline:
                        unique_map = 2 * int(eachline.strip().split()[0])
                    if 'aligned concordantly >1 times' in eachline:
                        mult_map = 2 * int(eachline.strip().split()[0])
                    if 'aligned discordantly 1 time' in eachline:
                        unique_map += 2 * int(eachline.strip().split()[0])
                    if 'aligned exactly 1 time' in eachline:
                        unique_map += int(eachline.strip().split()[0])
                    if 'aligned >1 times' in eachline:
                        mult_map += int(eachline.strip().split()[0])
                    if 'overall alignment rate' in eachline:
                        total_map_r = eachline.strip().split()[0]
                    if 'finished' in eachline:
                        sample = eachline.strip().split()[0]
                        total_map_n = '{0:.2f}'.format((unique_map + mult_map)/ 10**6)
                        mult_map_n = '{0:.2f}'.format(mult_map / 10**6)
                        mult_map_r = '{0:.2f}%'.format(100 * mult_map / total)
                        out_file_inf.write('{0}\t{1}\t{2}\n'.format(sample, total_map_r, mult_map_r))
        out_file_inf.close()





if __name__ == '__main__':
    main()
