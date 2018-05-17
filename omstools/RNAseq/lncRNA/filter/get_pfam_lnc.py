import pandas as pd
import click
import glob


def read_split_pfam(split_pfam_dir):
    pfam_df_list = list()
    split_pfam_out = glob.glob('{d}/*pfamA'.format(d=split_pfam_dir))
    for each_pfam in split_pfam_out:
        each_pfam_df = pd.read_table(each_pfam, header=None,
                                     delim_whitespace=True,
                                     skip_blank_lines=True,
                                     comment="#")
        pfam_df_list.append(each_pfam_df)
    pfam_df = pd.concat(pfam_df_list)
    return pfam_df


@click.command()
@click.argument(
    'split_pfam_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True,
)
@click.argument(
    'pfam_tucp',
    type=click.Path(dir_okay=False),
    required=True
)
@click.option(
    '-p',
    '--pfam_id',
    type=click.Path(dir_okay=False),
    help='valid pfam id list.',
    default=''
)
def main(split_pfam_dir, pfam_id, pfam_tucp):
    pfam_df = read_split_pfam(split_pfam_dir)
    pfam_df = pfam_df.loc[:, [0, 5]]
    pfam_df.columns = ['tr', 'pfam']
    if pfam_id:
        valid_pfam_df = pd.read_table(pfam_id, header=None, index_col=0)
        pfam_df = pfam_df[pfam_df.pfam.isin(valid_pfam_df.index)]
    pfam_df.loc[:, 'tr_id'] = [each.split('_')[0] for each in pfam_df.tr]
    pfam_df = pfam_df.drop_duplicates(subset=['tr_id'])
    pfam_df.to_csv(pfam_tucp, header=False,
                   index=False, columns=['tr_id'])


if __name__ == '__main__':
    main()
