import pandas as pd
from collections import Counter
import click
from scipy import stats
import glob


def read_split_pfam(split_pfam_dir):
    pfam_df_list = list()
    split_pfam_out = glob.glob('{d}/*pfamA'.format(d=split_pfam_dir))
    for each_pfam in split_pfam_out:
        try:
            each_pfam_df = pd.read_table(each_pfam, header=None,
                                         delim_whitespace=True,
                                         skip_blank_lines=True,
                                         comment="#")
        except pd.errors.EmptyDataError:
            continue
        else:
            pfam_df_list.append(each_pfam_df)
    pfam_df = pd.concat(pfam_df_list)
    return pfam_df


@click.command()
@click.option(
    '-c',
    '--coding_pfam',
    type=click.Path(exists=True),
    required=True,
    help='pfam prediction for coding sequences.')
@click.option(
    '-n',
    '--noncoding_pfam',
    type=click.Path(exists=True),
    required=True,
    help='pfam prediction for noncoding sequences.')
@click.option(
    '-p',
    '--pfam_hit',
    type=click.File('w'),
    required=True,
    help='valid Pfam domains in transcribed regions')
def main(coding_pfam, noncoding_pfam, pfam_hit):
    cd_df = read_split_pfam(coding_pfam)
    cd_total = len(cd_df)
    cd_count = Counter(cd_df.loc[:, 5])
    nc_df = read_split_pfam(noncoding_pfam)
    nc_total = len(nc_df)
    nc_count = Counter(nc_df.loc[:, 5])
    for each_hit in cd_count:
        if each_hit in nc_count:
            cd_hit, cd_nohit = cd_count[
                each_hit], cd_total - cd_count[each_hit]
            nc_hit, nc_nohit = nc_count[
                each_hit], nc_total - nc_count[each_hit]
            oddsratio, pvalue = stats.fisher_exact([[cd_hit, cd_nohit],
                                                    [nc_hit, nc_nohit]])
            if oddsratio > 10 and pvalue < 0.05:
                pfam_hit.write('{h}\n'.format(h=each_hit))
        else:
            pfam_hit.write('{h}\n'.format(h=each_hit))


if __name__ == '__main__':
    main()
