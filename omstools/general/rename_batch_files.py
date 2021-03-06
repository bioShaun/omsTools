import os
import pandas as pd
import click
import sys
from omstools.utils.config import CLICK_CONTEXT_SETTINGS
from omstools.utils.config import MutuallyExclusiveOption
reload(sys)
sys.setdefaultencoding('utf-8')


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.option(
    '-n',
    '--name_map',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
    help='Old name (first column) vs new name (second column).')
@click.option(
    '-s',
    '--suffix',
    type=click.STRING,
    default='',
    help='Suffix of files to rename.')
@click.option(
    '-p',
    '--prefix',
    type=click.STRING,
    default='',
    help='Prefix of files to rename.')
@click.option(
    '-d',
    '--file_dir',
    type=click.Path(file_okay=False, exists=True),
    required=True,
    help='Files directory.')
@click.option(
    '-i',
    '--inplace',
    is_flag=True,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["out_dir"],
    help='Rename files inplace.')
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['inplace'],
    help='Out directory to store renamed symbolic links of original files.')
@click.option(
    '--new_suffix',
    type=click.STRING,
    default=None,
    help='replace old suffix with new suffix.'
)
@click.option(
    '--new_prefix',
    type=click.STRING,
    default=None,
    help='replace old prefix with new prefix.'
)
def main(name_map, suffix, file_dir, inplace,
         out_dir, prefix, new_suffix, new_prefix):
    suf_len = len(suffix)
    pre_len = len(prefix)
    name_map_df = pd.read_table(name_map, delim_whitespace=True,
                                header=None, index_col=0)
    file_dir = os.path.abspath(file_dir)
    file_list = os.listdir(file_dir)
    for each_file in file_list:
        if prefix in each_file and suffix in each_file:
            each_name = each_file
            if prefix:
                each_name = each_name[pre_len:]
            if suffix:
                each_name = each_name[:-(suf_len)]
            each_file_path = os.path.join(file_dir, each_file)
            if each_name in name_map_df.index:
                new_name = name_map_df.loc[each_name, 1]
                out_suffix = suffix
                out_prefix = prefix
                if new_suffix:
                    out_suffix = new_suffix
                if new_prefix:
                    out_prefix = new_prefix
                if inplace:
                    new_file = os.path.join(file_dir, '{p}{n}{s}'.format(
                        n=new_name, s=out_suffix, p=out_prefix))
                    os.system('mv {o} {n}'.format(
                        o=each_file_path, n=new_file))
                else:
                    if not out_dir:
                        sys.exit(
                            'If not --inplace, --out_dir option is required.')
                    out_dir = os.path.abspath(out_dir)
                    new_file = os.path.join(out_dir, '{p}{n}{s}'.format(
                        n=new_name, s=out_suffix, p=out_prefix))
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    if not os.path.exists(new_file):
                        os.system(
                            'ln -s {o} {n}'.format(
                                o=each_file_path, n=new_file))
                        # print(
                        #     'ln -s {o} {n}'.format(
                        #         o=each_file_path, n=new_file))



if __name__ == '__main__':
    main()
