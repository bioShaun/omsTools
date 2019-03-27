import click
import os
import sys
from omstools.utils.systools import save_mkdir
from omstools.utils.config import gtf_tools
from omstools.utils.config import CLICK_CONTEXT_SETTINGS


def gtf_classifier(line, method):
    if 'attrs' in method:
        tag = method.split(':')[-1]
        if tag not in line.attrs:
            return 'unknown'
        else:
            tag_value = line.attrs[tag]
            if 'type' in tag:
                tag_value = gtf_tools['func_get_tr_type'](tag_value)
            return tag_value
    elif method == 'chr':
        return line.seqid
    else:
        sys.exit('unknown method')


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.option(
    '-g',
    '--gtf',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help='gtf file path.'
)
@click.option(
    '-m',
    '--method',
    type=click.STRING,
    default='attrs:gene_biotype',
    help='method to split gtf file, default is using \
    gene_biotype.supported method includes "chr" \
    and attribudes in your gtf file.'
)
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    default=os.getcwd(),
    help='directory to put splitted gtf files.'
)
def main(gtf, method, out_dir):
    '''Split gtf file according to its field or attribute
    '''
    save_mkdir(out_dir)
    gtf_split_dict = dict()
    with open(gtf) as gtf_inf:
        for eachline in gtf_tools['class_GTFFeature'].parse(gtf_inf):
            cat = gtf_classifier(eachline, method)
            gtf_split_dict.setdefault(cat, []).append(eachline)
    for each_cat in gtf_split_dict:
        each_cat_out = os.path.join(out_dir, '{c}.gtf'.format(c=each_cat))
        with open(each_cat_out, 'w') as gtf_out:
            for eachline in gtf_split_dict[each_cat]:
                gtf_out.write('{l_out}\n'.format(l_out=str(eachline)))


if __name__ == '__main__':
    main()
