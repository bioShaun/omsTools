import click
import os
from uti



@click.command()
@click.option(
    '-g',
    '--gtf',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    'm',
    '--method',
    type=click.STRING,
    default='attr:gene_biotype',
)
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd()
)
def main(gtf, method, out_dir):
    pass


if __name__ == '__main__':
    main()
