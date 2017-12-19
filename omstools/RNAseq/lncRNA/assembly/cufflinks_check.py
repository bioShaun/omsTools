import click
import pandas as pd
import glob
import os


def read_cufflinks_log(log_file):
    log_inf = open(log_file).readlines()
    for n, eachline in enumerate(log_inf):
        if 'Assembling transcripts and estimating abundances' in eachline:
            if (log_inf[n + 1]).strip()[-4:] == '100%':
                return True
    return False


@click.command()
@click.option(
    '-d',
    '--script_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True,
    help='cuflinks script directory.')
@click.option(
    '-o',
    '--output',
    type=click.File('w'),
    required=True,
    help='cufflinks run status.')
def main(script_dir, output):
    script_logs = glob.glob('{d}/*log.*'.format(d=script_dir))
    for each_log in script_logs:
        script_name = os.path.basename(each_log[:-19])
        if read_cufflinks_log(each_log):
            output.write('{n}\n'.format(n=script_name))


if __name__ == '__main__':
    main()
