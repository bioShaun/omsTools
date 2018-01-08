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
    required=True,
    help='output bed file path.'
)
@click.option(
    '-f',
    '--bed_format',
    default='bed12',
    type=click.Choice(['bed12', 'bed6', 'splicing', 'donor']),
    help='output bed format.'
)
def main(gtf, bed, bed_format):
    '''Convert gtf to bed12, bed6(including tss, exon and intron feature),
    splicing site and donor/acceptor site.
    '''
    gtf_dir_prefix = os.path.splitext(gtf.name)[0]
    try:
        transcript_objs = gtf_tools['func_tr_from_gtf_lines'](gtf)
    except GTFError:
        gtf.seek(0)
        formated_gtf_file = '{pfx}_formatted.gtf'.format(
            pfx=gtf_dir_prefix)
        gtf_tools['func_to_formatted_gtf'](gtf, formated_gtf_file)
        with open(formated_gtf_file) as gtf_inf:
            transcript_objs = gtf_tools['func_tr_from_gtf_lines'](gtf_inf)
    if bed_format == 'bed12':
        for each_tr_obj in transcript_objs:
            bed.write('{bedline}\n'.format(
                bedline=each_tr_obj.to_bed12()))
    elif bed_format == 'bed6':
        for each_tr_obj in transcript_objs:
            for each_line in each_tr_obj.to_feature_bed6():
                bed.write('{bedline}\n'.format(
                    bedline=each_line))
    elif bed_format == 'splicing':
        for each_tr_obj in transcript_objs:
            for eachline in each_tr_obj.to_splicing_region():
                bed.write('{bedline}\n'.format(
                    bedline=eachline))
    elif bed_format == 'donor':
        for each_tr_obj in transcript_objs:
            for eachline in each_tr_obj.to_donor_acceptor():
                bed.write('{bedline}\n'.format(
                    bedline=eachline))


if __name__ == '__main__':
    main()
