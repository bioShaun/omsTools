import click
import os
import pandas as pd
from omstools.utils.config import gtf_tools
from omstools.utils.config import CLICK_CONTEXT_SETTINGS
from omstools.utils.systools import save_mkdir


CURRENT_DIR = os.getcwd()
TR_HEADER = [
    'Transcript_id',
    'Gene_id',
    'Chr',
    'Start',
    'End',
    'Strand',
    'Exon_number',
    'Transcript_length'
]

EXON_HEADER = [
    'Chr',
    'Start',
    'End',
    'Strand',
    'Feature',
    'Transcript_id',
    'Gene_id',
    'Length'
]


def get_gene_feature(tr_df):
    gene_chr = tr_df.groupby(['Gene_id'])['Chr'].first()
    gene_start = tr_df.groupby(['Gene_id'])['Start'].min()
    gene_end = tr_df.groupby(['Gene_id'])['End'].max()
    gene_strand = tr_df.groupby(['Gene_id'])['Strand'].first()
    gene_tr = tr_df.groupby(['Gene_id'])['Transcript_id'].unique()
    gene_tr = pd.Series([','.join(each) for each in gene_tr],
                        index=gene_tr.index)
    gene_tr.name = 'Transcripts'
    gene_tr_num = tr_df.groupby(['Gene_id'])['Transcript_id'].count()
    gene_tr_num.name = 'Transcript_number'
    gene_exon = tr_df.groupby(['Gene_id'])['Exon_number'].max()
    gene_len = tr_df.groupby(['Gene_id'])['Transcript_length'].max()
    gene_df = pd.concat([gene_chr, gene_start, gene_end, gene_strand,
                        gene_tr, gene_tr_num, gene_exon, gene_len],
                        axis=1)
    gene_df = gene_df.reset_index()
    return gene_df


def add_type(func):
    def wrapper(gtf, biotype=False):
        if not biotype:
            return func(gtf)
        else:
            gene_feature_df, tr_feature_df, exon_intron_df = func(gtf)
            tr_type_dict = dict()
            for gene, tr_objs in gtf_tools['func_parse_gtf'](gtf):
                for each_tr in tr_objs:
                    tr_id = each_tr.attrs["transcript_id"]
                    tr_type = each_tr.attrs["transcript_biotype"]
                    gene_type = each_tr.attrs["gene_biotype"]
                    tr_type_dict.setdefault('Transcript_id', []).append(tr_id)
                    tr_type_dict.setdefault('Gene_id', []).append(gene)
                    tr_type_dict.setdefault(
                        'Transcript_biotype', []).append(tr_type)
                    tr_type_dict.setdefault(
                        'Gene_biotype', []).append(gene_type)
            tr_type_df = pd.DataFrame(tr_type_dict)
            gene_type_df = tr_type_df.loc[:, ['Gene_id', 'Gene_biotype']]
            tr_type_df = tr_type_df.drop('Gene_id', axis=1)
            gene_feature_df = pd.merge(
                gene_feature_df, gene_type_df,
                left_on='Gene_id', right_on='Gene_id', how='left')
            tr_feature_df = pd.merge(
                tr_feature_df, tr_type_df,
                left_on='Transcript_id', right_on='Transcript_id', how='left')
            exon_intron_df = pd.merge(
                exon_intron_df, tr_type_df,
                left_on='Transcript_id', right_on='Transcript_id', how='left')
            return gene_feature_df, tr_feature_df, exon_intron_df
    return wrapper


@add_type
def gtf2feature(gtf, biotype=False):
    tr_feature_dict = dict()
    exon_intron_dict = dict()
    for gene, tr_objs in gtf_tools['func_parse_gtf'](gtf):
        for each_tr in tr_objs:
            tr_id = each_tr.attrs["transcript_id"]
            exon_num = len(each_tr.exons)
            strand = gtf_tools['func_strand_int_to_str'](each_tr.strand)
            tr_len = each_tr.length
            tr_feature_dict.setdefault('Chr', []).append(each_tr.chrom)
            tr_feature_dict.setdefault('Start', []).append(each_tr.start + 1)
            tr_feature_dict.setdefault('End', []).append(each_tr.end)
            tr_feature_dict.setdefault('Strand', []).append(strand)
            tr_feature_dict.setdefault('Transcript_id', []).append(tr_id)
            tr_feature_dict.setdefault('Gene_id', []).append(gene)
            tr_feature_dict.setdefault('Exon_number', []).append(exon_num)
            tr_feature_dict.setdefault('Transcript_length', []).append(tr_len)
            for each_exon in each_tr.exons:
                each_exon_len = each_exon.end - each_exon.start
                exon_intron_dict.setdefault('Chr', []).append(each_tr.chrom)
                exon_intron_dict.setdefault('Start', []).append(
                    each_exon.start + 1)
                exon_intron_dict.setdefault('End', []).append(each_exon.end)
                exon_intron_dict.setdefault('Strand', []).append(strand)
                exon_intron_dict.setdefault('Feature', []).append('exon')
                exon_intron_dict.setdefault('Transcript_id', []).append(tr_id)
                exon_intron_dict.setdefault('Gene_id', []).append(gene)
                exon_intron_dict.setdefault('Length', []).append(each_exon_len)
            for each_intron in each_tr.iterintrons():
                intron_start, intron_end = each_intron
                intron_len = intron_end - intron_start
                intron_start += 1
                exon_intron_dict.setdefault('Chr', []).append(each_tr.chrom)
                exon_intron_dict.setdefault('Start', []).append(intron_start)
                exon_intron_dict.setdefault('End', []).append(intron_end)
                exon_intron_dict.setdefault('Strand', []).append(strand)
                exon_intron_dict.setdefault('Feature', []).append('intron')
                exon_intron_dict.setdefault('Transcript_id', []).append(tr_id)
                exon_intron_dict.setdefault('Gene_id', []).append(gene)
                exon_intron_dict.setdefault('Length', []).append(intron_len)
    gtf.seek(0)
    tr_feature_df = pd.DataFrame(tr_feature_dict)
    gene_feature_df = get_gene_feature(tr_feature_df)
    tr_feature_df = tr_feature_df.loc[:, TR_HEADER]
    exon_intron_df = pd.DataFrame(exon_intron_dict)
    exon_intron_df = exon_intron_df.loc[:, EXON_HEADER]
    return gene_feature_df, tr_feature_df, exon_intron_df


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.option(
    '-g',
    '--gtf',
    help='gtf file.',
    required=True,
    type=click.File('r')
)
@click.option(
    '-o',
    '--out_dir',
    help='directory to store assembly transcript/gene feature.',
    default=CURRENT_DIR,
    type=click.Path(file_okay=False)
)
@click.option(
    '--biotype',
    help='add transcript/gene biotype to feature table.',
    is_flag=True,
)
def main(gtf, out_dir, biotype):
    # insure out_dir exists
    save_mkdir(out_dir)
    # extract feature information from gtf file
    gene_feature_df, tr_feature_df, exon_intron_df = gtf2feature(
        gtf, biotype=biotype)
    # write out gene feature file
    gene_featrue_file = os.path.join(out_dir, 'Gene_feature.txt')
    gene_feature_df.to_csv(gene_featrue_file, sep='\t', index=False)
    # write out transcript feature file
    tr_feature_file = os.path.join(out_dir, 'Transcript_feature.txt')
    tr_feature_df.to_csv(tr_feature_file, sep='\t', index=False)
    # write out exon/intron feature file
    exon_intron_feature_file = os.path.join(out_dir, 'Exon_Intron_feature.txt')
    exon_intron_df.to_csv(exon_intron_feature_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
