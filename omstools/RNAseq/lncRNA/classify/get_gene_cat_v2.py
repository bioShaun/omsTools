import pandas as pd
import numpy as np
import click
import os


PRIORITY = ('Read-through', 'Protein coding',
            'Pseudogene', 'TUCP', 'lncrna', 'lncRNA', 'other', 'ncRNA,other')

type_map = {
    'other': 'lncRNA',
    'ncRNA,other': 'lncRNA',
    'lncrna': 'lncRNA',
    'protein_coding': 'Protein coding',
    'pseudogene': 'Pseudogene',
    'read_through': 'Read-through'
}


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
    '--out_dir',
    type=click.Path(file_okay=False),
    help='gene classify/summary directory based on \
    taco compare result and feelnc classify.',
    required=True
)
@click.option(
    '-n',
    '--name',
    type=click.STRING,
    help='Summary table name',
    default=None
)
def main(meta_table, tucp, out_dir, name):
    meta_table_df = pd.read_table(meta_table, index_col=0)
    tucp_df = pd.read_table(tucp, header=None, index_col=0)
    tucp_series = tucp_df.index.intersection(meta_table_df.index)
    meta_table_df.loc[tucp_series, 'category'] = 'TUCP'
    # label read_through
    for each_tr in meta_table_df.index:
        if meta_table_df.loc[each_tr, 'category_relative_detail'] == 'read_through':
            meta_table_df.loc[each_tr, 'category'] = 'read_through'
    meta_table_df = meta_table_df[meta_table_df.category_relative_detail !=
                                  'intronic_same_strand']
    meta_table_df.loc[:, 'category'].replace(type_map, inplace=True)
    meta_table_df = meta_table_df.reset_index()
    tr_type_summary = meta_table_df.category.value_counts()
    meta_table_df.loc[:, 'novel_status'] = np.where(
        meta_table_df.category_relative == 'exonic_overlap', 'known', 'novel')
    meta_table_df.loc[:, 'new_gene_id'] = meta_table_df.novel_status + \
        '.' + meta_table_df.gene_id
    meta_table_type_df = meta_table_df.loc[:, ['new_gene_id', 'category']]
    meta_table_type_df.columns = ['gene_id', 'category']
    gene_type_map = meta_table_type_df.groupby(
        ['gene_id'])['category'].unique()
    meta_table_df = meta_table_df.reset_index()
    gene_name_df = meta_table_df.loc[:, ['new_gene_id',
                                         'category_relative',
                                         'ref_gene_id',
                                         'ref_gene_name']]
    gene_name_df.columns = [
        'gene_id', 'category_relative', 'ref_gene_id', 'ref_gene_name']
    gene_name_df = gene_name_df[gene_name_df.category_relative ==
                                'exonic_overlap']
    gene_name_df = gene_name_df.loc[:, [
        'gene_id', 'ref_gene_id', 'ref_gene_name']].drop_duplicates()

    def get_type(type_list):
        for each_type in PRIORITY:
            if each_type in type_list:
                return type_map.get(each_type, each_type)

    gene_type_list = map(get_type, gene_type_map)
    gene_type_df = pd.DataFrame(
        gene_type_list, index=gene_type_map.index, columns=['type'])
    read_through_genes = gene_type_df[gene_type_df.type ==
                                      "read_through"].index
    gene_name_df = gene_name_df[~gene_name_df.gene_id.isin(read_through_genes)]
    gene_name_df = gene_name_df.set_index('gene_id')
    read_through_sup = gene_name_df[
        gene_name_df.index.value_counts() > 1].index.unique()
    gene_type_df.loc[read_through_sup, 'type'] = 'Read-through'
    gene_type_summary = gene_type_df.type.value_counts()
    type_stats = pd.concat([tr_type_summary, gene_type_summary], axis=1)
    type_stats.columns = ['Transcript', 'Gene']
    type_stats.index.name = 'Category'
    summary_file = os.path.join(out_dir, 'assembly.number.summary.txt')
    classify_file = os.path.join(out_dir, 'gene.classify.txt')
    name_file = os.path.join(out_dir, 'gene.name.txt')
    if name is not None:
        type_stats.loc[:, 'Name'] = name
    type_stats.to_csv(summary_file, sep='\t', header=False)
    gene_type_df.to_csv(classify_file, sep='\t')
    gene_name_df = gene_name_df[gene_name_df.index.value_counts() == 1]
    gene_name_df.to_csv(name_file, sep='\t')


if __name__ == '__main__':
    main()
