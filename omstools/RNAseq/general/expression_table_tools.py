#! /usr/bin/env python

import pandas as pd
import click
'''
gene expression matrix, with gene id in first column,
gene expression level of each sample in othre columns.
'''


@click.group(chain=True, invoke_without_command=True)
@click.argument('exp_table', type=click.STRING, required=True)
@click.pass_context
def main(ctx, exp_table):
    ctx.obj['exp_table'] = exp_table


@main.command('merge_by_group')
@click.option(
    '-s',
    '--sample_inf',
    type=click.STRING,
    required=True,
    help='sample vs group file, with group id in first column,\
              sample id in second column, seperated with tab.')
@click.option(
    '-o',
    '--output',
    type=click.STRING,
    default='genes.group.matrix.txt',
    help='table with mean expression level of each group.')
@click.pass_context
def merge_by_group(ctx, sample_inf, output):
    sample_df = pd.read_table(sample_inf, header=None, index_col=1)
    gene_exp_df = pd.read_table(ctx.obj['exp_table'], index_col=0)
    sample_df.columns = ['Group']
    merged_df = pd.merge(
        sample_df, gene_exp_df.T, left_index=True, right_index=True)
    merged_df_group = merged_df.groupby(['Group'])
    out_df = merged_df_group.mean().T
    out_df.to_csv(output, sep='\t')


if __name__ == '__main__':
    main(obj={})
