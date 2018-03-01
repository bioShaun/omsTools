import pandas as pd
import click


PRIORITY = ('mixed_read_through', 'protein_coding',
            'pseudogene', 'TUCP', 'lncrna', 'lncRNA', 'other', 'ncRNA,other')


@click.command()
@click.option(
    '-m',
    '--meta_table',
    type=click.Path(exists=True, dir_okay=False),
    help='taco compare metadata',
    required=True,
)
@click.option(
    '-t',
    '--tucp',
    type=click.Path(exists=True, dir_okay=False),
    help='tucp transcripts.',
    required=True,
)
@click.option(
    '-o',
    '--output',
    type=click.Path(dir_okay=False),
    help='gene classify based on taco compare result and feelnc classify.',
    default='gene.classify.txt'
)
def main(meta_table, tucp, output):
    meta_table_df = pd.read_table(meta_table, index_col=0)
    meta_table_type_df = meta_table_df.loc[:, ['gene_id', 'category']]
    tucp_df = pd.read_table(tucp, header=None, index_col=0)
    tucp_detail = meta_table_df.loc[tucp_df.index]
    tucp_detail = tucp_detail.loc[:, ['gene_id', 'category']]
    tucp_detail.loc[:, 'category'] = 'TUCP'
    merged_type_table = pd.concat([meta_table_type_df, tucp_detail])
    gene_type_map = merged_type_table.groupby(['gene_id'])['category'].unique()

    def get_type(type_list):
        for each_type in PRIORITY:
            if each_type in type_list:
                if 'other' in each_type:
                    each_type = 'lncrna'
                return each_type

    gene_type_list = map(get_type, gene_type_map)
    gene_type_df = pd.DataFrame(
        gene_type_list, index=gene_type_map.index, columns=['type'])
    gene_type_df.to_csv(output, sep='\t')


if __name__ == '__main__':
    main()
