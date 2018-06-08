import pandas as pd
from collections import Counter
import click
from scipy import stats
import glob
import os
import math


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
    '-o',
    '--pfam_out',
    type=click.Path(file_okay=False),
    required=True,
    help='Pfam analysis output directory.')
def main(coding_pfam, noncoding_pfam, pfam_out):
    cd_df = read_split_pfam(coding_pfam)
    cd_total = len(cd_df)
    cd_pfam = len(cd_df.loc[:, 5].unique())
    cd_count = Counter(cd_df.loc[:, 5])
    nc_df = read_split_pfam(noncoding_pfam)
    nc_total = len(nc_df)
    nc_pfam = len(nc_df.loc[:, 5].unique())
    nc_count = Counter(nc_df.loc[:, 5])
    valid_pfam_list = list()
    artif_pfam_list = list()
    pfam_hit_file = os.path.join(pfam_out, 'valid.pfam')
    pfam_hit = open(pfam_hit_file, 'w')
    pfam_plot_file = os.path.join(pfam_out, 'pfam.detail.txt')
    pfam_plot = open(pfam_plot_file, 'w')
    pfam_plot.write(
        "pfam\tcoding_region_hit\tnon-transcribed_region_hit\tother_domain_coding_region_hit\tother_domain_non-transcribed_region_hit\todds_ratio\tpvalue\t\tstatus\n")
    for each_hit in cd_count:
        cd_hit, cd_nohit = cd_count[
            each_hit], cd_total - cd_count[each_hit]
        is_valid = 'Valid Pfam domains'
        if each_hit in nc_count:
            nc_hit, nc_nohit = nc_count[
                each_hit], nc_total - nc_count[each_hit]
        else:
            nc_hit = 0
            nc_nohit = nc_total
        oddsratio, pvalue = stats.fisher_exact([[cd_hit, nc_hit],
                                                [cd_nohit, nc_nohit]])
        # if math.isinf(oddsratio):
        #     pvalue = 0
        if oddsratio > 10 and pvalue < 0.05:
            valid_pfam_list.append(each_hit)
            pfam_hit.write('{h}\n'.format(h=each_hit))
        else:
            artif_pfam_list.append(each_hit)
            is_valid = 'Artifact Pfam domains'
        pfam_plot.write('{hit}\t{ch}\t{nh}\t{cn}\t{nn}\t\t{od}\t{pv}\t{vl}\n'.format(
            hit=each_hit, ch=cd_hit,
            nh=nc_hit, vl=is_valid,
            cn=cd_nohit, nn=nc_nohit,
            od=oddsratio, pv=pvalue
        ))
    pfam_hit.close()
    pfam_plot.close()
    valid_pfam_num = len(valid_pfam_list)
    artif_pfam_num = len(artif_pfam_list)
    valid_pfam_hits = len(cd_df[cd_df.loc[:, 5].isin(valid_pfam_list)])
    artif_pfam_hits = cd_total - valid_pfam_hits
    nc_valid_hits = len(nc_df[nc_df.loc[:, 5].isin(valid_pfam_list)])
    print '''
Coding region pfam hits:             {ch}
Coding region pfam domains:          {cd}
Non-transcribed region pfam hits:    {nh}
Non-transcribed region pfam domains: {nd}
Valid pfam domains:                  {vd}
Artifacts pfam domains:              {ad}
Valid pfam hits:                     {vh}
Artifacts pfam hits:                 {ah}
Non-transcribed region valid pfam hits: {nvh}
'''.format(ch=cd_total, cd=cd_pfam, nh=nc_total, nd=nc_pfam,
           vd=valid_pfam_num, ad=artif_pfam_num, vh=valid_pfam_hits,
           ah=artif_pfam_hits, nvh=nc_valid_hits)


if __name__ == '__main__':
    main()
