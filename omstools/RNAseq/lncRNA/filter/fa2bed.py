import click
from Bio import SeqIO


@click.command()
@click.argument(
    'fasta',
    type=click.Path(dir_okay=False, exists=True),
    required=True
)
@click.argument(
    'bed',
    type=click.File('w'),
    required=True
)
def main(fasta, bed):
    chrom = '1'
    for seq_record in SeqIO.parse(fasta, "fasta"):
        bed.write('{chrom}\t0\t{end}\t{s_id}\t.\t+\n'.format(
            chrom=chrom, end=len(seq_record.seq),
            s_id=seq_record.id
        ))


if __name__ == '__main__':
    main()
