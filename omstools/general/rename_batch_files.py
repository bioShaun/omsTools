import os
import pandas as pd
import click
import sys
reload(sys)
sys.setdefaultencoding('utf-8')


class MutuallyExclusiveOption(click.Option):
    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop('mutually_exclusive', []))
        help = kwargs.get('help', '')
        if self.mutually_exclusive:
            ex_str = ', '.join(self.mutually_exclusive)
            kwargs['help'] = help + (
                ' NOTE: This argument is mutually exclusive with '
                ' arguments: [' + ex_str + '].'
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise click.UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(
                    self.name,
                    ', '.join(self.mutually_exclusive)
                )
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(
            ctx,
            opts,
            args
        )


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
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
def main(name_map, suffix, file_dir, inplace, out_dir, prefix):
    suf_len = len(suffix)
    pre_len = len(prefix)
    name_map_df = pd.read_table(name_map, header=None, index_col=0)
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
                if inplace:
                    new_file = os.path.join(file_dir, '{p}{n}{s}'.format(
                        n=new_name, s=suffix, p=prefix))
                    os.system('mv {o} {n}'.format(
                        o=each_file_path, n=new_file))
                else:
                    if not out_dir:
                        sys.exit(
                            'If not --inplace, --out_dir option is required.')
                    out_dir = os.path.abspath(out_dir)
                    new_file = os.path.join(out_dir, '{p}{n}{s}'.format(
                        n=new_name, s=suffix, p=prefix))
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    if not os.path.exists(new_file):
                        os.system(
                            'ln -s {o} {n}'.format(
                                o=each_file_path, n=new_file))
                    # print 'ln -s {o} {n}'.format(o=each_file, n=new_file)


if __name__ == '__main__':
    main()
