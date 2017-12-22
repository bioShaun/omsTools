import click
import os
from omstools.utils.config import gtf_tools
from omstools.utils.config import CLICK_CONTEXT_SETTINGS

GTFError = gtf_tools['error_GTFError']


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.option(
    '-g',
    '--gtf',
    type=click.File('r'),
    required=True,
    help='gtf file path.'
)
@click.option(
    '-b',
    '--bed',
    type=click.File('w'),
    default='',
    help='output bed file path.'
)
@click.option(
    '-f',
    '--format',
    default='bed12',
    type=click.Choice(['bed12', 'bed6']),
    help='output bed format.'
)
def main(gtf, bed, format):
    gtf_dir_prefix = os.path.splitext(gtf.name)[0]
    try:
        transcript_objs = gtf_tools['func_tr_from_gtf_lines'](gtf)
    except GTFError:
        gtf.seek(0)
        formated_gtf_file = '{pfx}_formatted.gtf'.format(
            pfx=gtf_dir_prefix)
        gtf_tools['func_to_formatted_gtf'](gtf, formated_gtf_file)


if __name__ == '__main__':
    main()
