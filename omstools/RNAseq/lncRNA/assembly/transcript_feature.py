import click
from transcript import parse_gtf, strand_int_to_str
import os
import pandas as pd


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
    'Gene_id',
    'Chr',
    'Start',
    'End',
    'Strand',
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
    gene_tr = pd.Series([','.join(each) for each in gene_chr],
                        index=gene_chr.index)
    gene_tr.name = 'Transcripts'
    gene_tr_num = tr_df.groupby(['Gene_id'])['Transcript_id'].count()
    gene_tr_num.name = 'Transcript_number'
    gene_exon = tr_df.groupby(['Gene_id'])['Exon_number'].max()
    gene_len = tr_df.groupby(['Gene_id'])['Transcript_length'].max()
    gene_df = pd.concat([gene_chr, gene_start, gene_end, gene_strand,
                        gene_tr, gene_tr_num, gene_exon, gene_len],
                        axis=1)
    return gene_df


def gtf2feature(gtf):
    tr_feature_dict = dict()
    exon_intron_dict = dict()
    gtf_inf = open(gtf, 'r')
    for gene, tr_objs in parse_gtf(gtf_inf):
        for each_tr in tr_objs:
            tr_id = each_tr.attrs["transcript_id"]
            exon_num = len(each_tr.exons)
            strand = strand_int_to_str(each_tr.strand)
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
                exon_intron_dict.setdefault('Type', []).append('exon')
                exon_intron_dict.setdefault('Transcript_id', []).append(tr_id)
                exon_intron_dict.setdefault('Gene_id', []).append(gene)
                exon_intron_dict.setdefault('Length', []).append(each_exon_len)
    gtf_inf.close()
    tr_feature_df = pd.DataFrame(tr_feature_dict)
    gene_feature_df = get_gene_feature(tr_feature_df)
    tr_feature_df = tr_feature_df.loc[:, TR_HEADER]
    exon_intron_df = pd.DataFrame(exon_intron_dict)
    exon_intron_df = exon_intron_df.loc[:, EXON_HEADER]
    return gene_feature_df, tr_feature_df, exon_intron_df


@click.command()
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
    type=click.File('w')
)
def main(gtf, out_dir):
    pass


if __name__ == '__main__':
    main()
