import click
import os
import pandas as pd
from omstools.utils.config import gtf_tools
from omstools.utils.config import CLICK_CONTEXT_SETTINGS
from omstools.utils.systools import save_mkdir
import envoy

script_dir, script_name = os.path.split(os.path.abspath(__file__))
FEATURE_PLOT_R = os.path.join(script_dir, 'transcript_feature.R')

CURRENT_DIR = os.getcwd()
TR_HEADER = [
    'Transcript_id', 'Gene_id', 'Chr', 'Start', 'End', 'Strand', 'Exon_number',
    'Transcript_length'
]

EXON_HEADER = [
    'Chr', 'Start', 'End', 'Strand', 'Feature', 'Transcript_id', 'Gene_id',
    'Length'
]


def is_formatted(gtf):
    '''check if gtf file is formated
    '''
    gtf_obj = gtf_tools['func_parse_gtf'](gtf)
    try:
        gtf_obj.next()
    except gtf_tools['error_GTFError']:
        flag = False
    else:
        flag = True
    finally:
        gtf.seek(0)
        return flag


def format_gtf(gtf):
    '''format gtf file to assembly line accepted format
    '''
    if is_formatted(gtf):
        return gtf
    else:
        gtf_path, gtf_name = os.path.split(gtf.name)
        formated_gtf = os.path.join(
            gtf_path, 'formated.{n}'.format(n=gtf_name))
        gtf_tools['func_to_formatted_gtf'](gtf, formated_gtf)
        gtf_inf = open(formated_gtf, 'r')
        return gtf_inf


def get_gene_feature(tr_df):
    '''summarize gene feature according to transcript feature table
    '''
    gene_chr = tr_df.groupby(['Gene_id'])['Chr'].first()
    gene_start = tr_df.groupby(['Gene_id'])['Start'].min()
    gene_end = tr_df.groupby(['Gene_id'])['End'].max()
    gene_strand = tr_df.groupby(['Gene_id'])['Strand'].first()
    gene_tr = tr_df.groupby(['Gene_id'])['Transcript_id'].unique()
    gene_tr = pd.Series(
        [','.join(each) for each in gene_tr], index=gene_tr.index)
    gene_tr.name = 'Transcripts'
    gene_tr_num = tr_df.groupby(['Gene_id'])['Transcript_id'].count()
    gene_tr_num.name = 'Transcript_number'
    gene_exon = tr_df.groupby(['Gene_id'])['Exon_number'].max()
    gene_len = tr_df.groupby(['Gene_id'])['Transcript_length'].max()
    gene_df = pd.concat(
        [
            gene_chr, gene_start, gene_end, gene_strand, gene_tr, gene_tr_num,
            gene_exon, gene_len
        ],
        axis=1)
    gene_df = gene_df.reset_index()
    return gene_df


def correct_key_val(mydict, key_list):
    val = None
    for each_key in key_list:
        if mydict.get(each_key) is not None:
            val = mydict.get(each_key)
            return val
    if val is None:
        raise KeyError


def add_type(func):
    '''add transcript/gene biotype to gene/transcript feature table
    '''

    def wrapper(gtf, biotype=False):
        if not biotype:
            return func(gtf)
        else:
            gene_feature_df, tr_feature_df, exon_intron_df = func(gtf)
            tr_type_dict = dict()
            for gene, tr_objs in gtf_tools['func_parse_gtf'](gtf):
                for each_tr in tr_objs:
                    tr_id = each_tr.attrs["transcript_id"]
                    tr_type = correct_key_val(
                        each_tr.attrs, ["transcript_biotype", 'transcript_type'])
                    gene_type = correct_key_val(
                        each_tr.attrs, ["gene_biotype", 'gene_type'])
                    gene_type = gtf_tools['func_get_tr_type'](gene_type)
                    if gene_type != 'lncRNA':
                        tr_type = gtf_tools['func_get_tr_type'](tr_type)
                    tr_type_dict.setdefault('Transcript_id', []).append(tr_id)
                    tr_type_dict.setdefault('Gene_id', []).append(gene)
                    tr_type_dict.setdefault('Transcript_biotype',
                                            []).append(tr_type)
                    tr_type_dict.setdefault('Gene_biotype',
                                            []).append(gene_type)
            tr_type_df = pd.DataFrame(tr_type_dict)
            gene_type_df = tr_type_df.loc[:, ['Gene_id', 'Gene_biotype']]
            gene_type_df = gene_type_df.drop_duplicates()
            tr_type_df = tr_type_df.drop('Gene_id', axis=1)
            gene_feature_df = pd.merge(
                gene_feature_df,
                gene_type_df,
                left_on='Gene_id',
                right_on='Gene_id',
                how='left')
            tr_feature_df = pd.merge(
                tr_feature_df,
                tr_type_df,
                left_on='Transcript_id',
                right_on='Transcript_id',
                how='left')
            exon_intron_df = pd.merge(
                exon_intron_df,
                tr_type_df,
                left_on='Transcript_id',
                right_on='Transcript_id',
                how='left')
            return gene_feature_df, tr_feature_df, exon_intron_df

    return wrapper


# def sort_table(func):
#     def wrapper(gtf):
#         pass

#     pass

@add_type
def gtf2feature(gtf, biotype=False):
    '''extract transcript/exon/intron basic information from gtf
    and store it to DataFrame
    '''
    tr_feature_dict = dict()
    exon_intron_dict = dict()
    gtf_obj_iter = gtf_tools['func_parse_gtf'](gtf)
    for gene, tr_objs in gtf_obj_iter:
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
                exon_intron_dict.setdefault('Start',
                                            []).append(each_exon.start + 1)
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


def plot_feature(out_dir):
    '''generate feature plot in output directory
    '''
    cmd1 = 'Rscript {rs} --feature_file_dir {dr} --out_dir {dr}'.format(
        rs=FEATURE_PLOT_R, dr=out_dir)
    cmd2 = '{cmd} --detail'.format(cmd=cmd1)
    cmd_list = [cmd1, cmd2]
    err_list = []
    for each_cmd in cmd_list:
        r = envoy.run(each_cmd)
        print r.std_err
        err_list.append(r.std_err)
    return err_list


def get_summary(df, cat_col, stat_col):
    summary_df_list = list()
    for each_stat in stat_col:
        each_stat_list = list()
        df_group = df.groupby(cat_col)[each_stat]
        for i in (.25, .5, .75):
            each_stat_list.append(df_group.quantile(i))
        each_stat_list.append(df_group.mean())
        each_stat_df = pd.concat(each_stat_list, axis=1)
        each_stat_df.loc[:, 'Stat'] = each_stat
        each_stat_df.columns = ['Q1', 'Median', 'Q3', 'Mean', 'Stat']
        summary_df_list.append(each_stat_df)
    summary_df = pd.concat(summary_df_list)
    return summary_df


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.option(
    '-g',
    '--gtf',
    help='gtf file.',
    required=True,
    type=click.File('r'))
@click.option(
    '-o',
    '--out_dir',
    help='directory to store assembly transcript/gene feature.',
    default=CURRENT_DIR,
    type=click.Path(file_okay=False))
@click.option(
    '--biotype',
    help='add transcript/gene biotype to feature table.',
    is_flag=True,
)
@click.option(
    '--plot',
    help='generate feature plot.',
    is_flag=True
)
def main(gtf, out_dir, biotype, plot):
    '''Extract basic information of gene/transcript/exon/intron from gtf file
    '''
    # insure out_dir exists
    save_mkdir(out_dir)
    # extract feature information from gtf file
    gene_featrue_file = os.path.join(out_dir, 'Gene_feature.txt')
    tr_feature_file = os.path.join(out_dir, 'Transcript_feature.txt')
    exon_intron_feature_file = os.path.join(out_dir, 'Exon_Intron_feature.txt')
    gtf = format_gtf(gtf)
    gene_feature_df, tr_feature_df, exon_intron_df = gtf2feature(
        gtf, biotype=biotype)

    # write out gene feature file
    gene_feature_df.to_csv(gene_featrue_file, sep='\t', index=False)

    # filter tr_feature_df lncRNA in protein coding gene
    tr_feature_df = tr_feature_df[~((tr_feature_df.Gene_biotype == "protein_coding") &
                                    (tr_feature_df.Transcript_biotype != "protein_coding"))]
    exon_intron_df = exon_intron_df[~((exon_intron_df.Gene_biotype == "protein_coding") &
                                      (exon_intron_df.Transcript_biotype != "protein_coding"))]
    # write out transcript feature file
    tr_feature_df.to_csv(tr_feature_file, sep='\t', index=False)
    # write out exon/intron feature file
    exon_intron_df.to_csv(exon_intron_feature_file, sep='\t', index=False)
    # plot output
    if plot:
        plot_feature(out_dir)
    # output summary
    if biotype:
        gene_feature_summary = os.path.join(
            out_dir, 'Gene_feature_summary.txt')
        tr_feature_summary = os.path.join(
            out_dir, 'Transcript_feature_summary.txt')
        exon_intron_feature_summary = os.path.join(
            out_dir, 'Exon_Intron_feature_summary.txt')
        gene_summary_df = get_summary(
            gene_feature_df, 'Gene_biotype',
            ['Transcript_number', 'Exon_number', 'Transcript_length'])
        gene_summary_df.to_csv(gene_feature_summary, sep='\t')

        # def tr_feature_summary(df, , stat_col):
        #     df.loc[:, 'Transcript_biotype_general'] = map(
        #         gtf_tools['func_get_tr_type'],
        #         tr_feature_df.loc[:, 'Transcript_biotype'])
        #     summary_df1 = get_summary(
        #         tr_feature_df, 'Transcript_biotype',
        #         stat_col)
        #     summary_df2 = get_summary(
        #         tr_feature_df, 'Transcript_biotype_general',
        #         stat_col)
        #     summary_df = pd.concat([summary_df1, summary_df2])
        #     summary_df = summary_df.drop_duplicates()

        tr_feature_df.loc[:, 'Transcript_biotype_general'] = map(
            gtf_tools['func_get_tr_type'],
            tr_feature_df.loc[:, 'Transcript_biotype'])
        tr_summary_df1 = get_summary(
            tr_feature_df, 'Transcript_biotype',
            ['Exon_number', 'Transcript_length'])
        tr_summary_df2 = get_summary(
            tr_feature_df, 'Transcript_biotype_general',
            ['Exon_number', 'Transcript_length'])
        tr_summary_df = pd.concat([tr_summary_df1, tr_summary_df2])
        tr_summary_df = tr_summary_df.drop_duplicates()
        tr_summary_df.index.name = 'Transcript_biotype'
        tr_summary_df = tr_summary_df.reset_index()
        tr_summary_df = tr_summary_df.sort_values(
            by=['Stat', 'Transcript_biotype'])
        tr_summary_df.to_csv(tr_feature_summary, sep='\t', index=False)
        ei_id_tuple = zip(exon_intron_df.Chr, exon_intron_df.Start.map(str),
                          exon_intron_df.End.map(str), exon_intron_df.Strand,
                          exon_intron_df.Transcript_biotype)
        ei_id_list = ['_'.join(each) for each in ei_id_tuple]
        exon_intron_df.loc[:, 'exon_id'] = ei_id_list
        exon_intron_df.loc[:, 'Transcript_biotype_general'] = map(
            gtf_tools['func_get_tr_type'],
            exon_intron_df.loc[:, 'Transcript_biotype'])
        ei_summary_df1 = get_summary(
            exon_intron_df,
            ['Transcript_biotype', 'Feature'], ['Length'])
        ei_summary_df2 = get_summary(
            exon_intron_df,
            ['Transcript_biotype_general', 'Feature'], ['Length'])
        ei_summary_df = pd.concat([ei_summary_df1, ei_summary_df2])
        ei_summary_df = ei_summary_df.drop_duplicates().sort_index()
        ei_summary_df.to_csv(exon_intron_feature_summary, sep='\t')
    # TODO sort output table


if __name__ == '__main__':
    main()
