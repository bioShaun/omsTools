#! /urs/bin/env python

from __future__ import division
import click
import os
from HTSeq import GFF_Reader
import pandas as pd
import sys


INTERGENIC_INF = ['-', '-', '-', '-', '-', '-', '-', 'intergenic']
LNC_ORDER = ['divergent', 'sense_intronic',
             'other_sense_overlap', 'antisense', 'intergenic']
OVERLAP_CUTOFF = 0.5


def iterintrons(tss, exon_sizes, exon_starts):
    tss = int(tss)
    exon_sizes_list = [int(each) for each in exon_sizes.split(',') if each]
    exon_starts_list = [int(each) for each in exon_starts.split(',') if each]
    exon_list = list()
    for n, each_size in enumerate(exon_sizes_list):
        each_start = exon_starts_list[n]
        exon_list.append((tss + each_start,
                          tss + each_start + each_size))
    e1 = exon_list[0]
    for n in range(1, len(exon_list)):
        e2 = exon_list[n]
        yield e1[1], e2[0]
        e1 = e2


def tss_in_interval(tss, iter_interval):
    for each in iter_interval:
        if tss >= each[0] and tss < each[1]:
            return True
    else:
        return False


def overlap_portion(inter_rd):
    overlap_len = inter_rd[24]
    tr1_len = inter_rd[2] - inter_rd[1]
    tr2_len = inter_rd[14] - inter_rd[13]
    return tr1_len / overlap_len, tr2_len / overlap_len


def compare_class(lnc_class1, lnc_class2):
    compare1 = lnc_class1[:]
    compare2 = lnc_class2[:]
    compare1[0] = list(reversed(LNC_ORDER)).index(compare1[0])
    compare2[0] = list(reversed(LNC_ORDER)).index(compare2[0])
    compare1[1] *= -1
    compare2[1] *= -1
    compare_list = [compare1, compare2]
    compare_list.sort()
    return compare_list.index(compare1)


@click.command()
@click.option(
    '-g',
    '--gtf',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help='lncRNA candidates gtf.')
@click.option(
    '-f',
    '--feelnc_classify',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help='FEElnc classify result.')
@click.option(
    '-b',
    '--bed_intersect',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help='bedtools intersect for lncRNA bed and mRNA bed.')
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd(),
    help='output directory.')
def main(gtf, feelnc_classify, bed_intersect, out_dir):
    feelnc_df = pd.read_table(feelnc_classify, index_col=2)
    intersect_df = pd.read_table(bed_intersect, index_col=[3, 15], header=None)
    lnc_class_list = []
    out_header = list(feelnc_df.columns[1:])
    out_header.insert(0, 'lncRNA_transcript')
    out_header.append('lncRNA_class')

    def get_class(fee_rd, intersect_df):
        if fee_rd.type == 'intergenic':
            if fee_rd.subtype == 'divergent':
                return 'divergent', 0, 0
            else:
                return 'intergenic', 0, 0
        else:
            inter_index = (fee_rd.name, fee_rd.partnerRNA_transcript)
            inter_rd = intersect_df.loc[inter_index]
            overlap1, overlap2 = overlap_portion(inter_rd)
            if fee_rd.direction == 'sense':
                if fee_rd.subtype == 'containing':
                    return 'other_sense_overlap', overlap1, overlap2
                elif fee_rd.subtype == 'nested':
                    return 'sense_intronic', overlap1, overlap2
                elif fee_rd.subtype == 'overlapping':
                    if overlap1 >= OVERLAP_CUTOFF:
                        introns = iterintrons(inter_rd[13], inter_rd[22],
                                              inter_rd[23])
                        lnc_can_start = inter_rd[1]
                        if tss_in_interval(lnc_can_start, introns):
                            return 'sense_intronic', overlap1, overlap2
                    return 'other_sense_overlap', overlap1, overlap2
                else:
                    sys.exit('unkown type [{t.subtype}]'.format(t=fee_rd))
            else:
                if fee_rd.subtype == 'nested':
                    return 'antisense', overlap1, overlap2
                else:
                    if overlap1 >= OVERLAP_CUTOFF:
                        return 'antisense', overlap1, overlap2
                    else:
                        return 'intergenic', overlap1, overlap2

    def lnc_classify(tr_id, feelnc_df, intersect_df):
        tr_detail_df = feelnc_df.loc[tr_id]
        out_inf = []
        class_inf = []
        if tr_detail_df.index[0] == tr_id:
            for n in range(len(tr_detail_df)):
                class_value = list(get_class(tr_detail_df.ix[n], intersect_df))
                dis = tr_detail_df.ix[n].distance
                tmp_class_inf = class_value[:]
                tmp_class_inf.insert(1, dis)
                tmp_out_inf = list(tr_detail_df.ix[n][1:])
                if not out_inf:
                    out_inf = tmp_out_inf
                    class_inf = tmp_class_inf
                else:
                    if compare_class(tmp_class_inf, class_inf):
                        out_inf = tmp_out_inf
                        class_inf = tmp_class_inf

        else:
            class_value = list(get_class(tr_detail_df, intersect_df))
            out_inf = list(tr_detail_df[1:])
            class_inf = class_value
        out_inf.insert(0, tr_id)
        out_inf.append(class_inf[0])
        return out_inf

    for eachline in GFF_Reader(gtf):
        if eachline.type == 'transcript':
            tr_id = eachline.attr['transcript_id']
            gene_id = eachline.attr['gene_id']
            if tr_id not in feelnc_df.index:
                out_inf = [tr_id, gene_id]
                out_inf.extend(INTERGENIC_INF)
            else:
                out_inf = lnc_classify(tr_id, feelnc_df, intersect_df)
            out_inf_series = pd.Series(out_inf, index=out_header)
            lnc_class_list.append(out_inf_series)

    out_df = pd.concat(lnc_class_list, axis=1).T
    out_file = os.path.join(out_dir, 'lncRNA.classify.txt')
    out_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
