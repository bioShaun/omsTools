import click
import gtfparse
from omstools.utils.config import gtf_tools


@click.command()
@click.argument(
    'gtf',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'coding_region_bed',
    type=click.Path(exists=False),
    required=True
)
def main(gtf, coding_region_bed):
    gtf_df = gtfparse.read_gtf(gtf)
    gene_df = gtf_df[gtf_df.feature == 'gene']
    gene_df.gene_biotype.replace(gtf_tools['dict_GENCODE_CATEGORY_MAP'],
                                 inplace=True)
    gene_df = gene_df[(gene_df.gene_biotype == 'protein_coding') |
                      (gene_df.gene_biotype == 'pseudogene')]
    coding_region_bed_df = gene_df.loc[:, [
        'seqname', 'start', 'end', 'gene_id', 'score', 'strand']]
    coding_region_bed_df.to_csv(coding_region_bed, header=None,
                                index=False, na_rep='.', sep='\t')


if __name__ == '__main__':
    main()
