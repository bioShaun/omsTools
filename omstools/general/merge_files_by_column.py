import click
import pandas as pd
import os


def merge_files(df_list, df0=None):
    if df0 is None:
        df0 = df_list.pop(0)
        return merge_files(df_list, df0)
    elif df_list:
        df1 = df_list.pop(0)
        new_df = pd.merge(df0, df1,
                          left_index=True,
                          right_index=True)
        new_df.index.name = df0.index.name
        return merge_files(df_list, df0=new_df)
    else:
        return df0


@click.command()
@click.argument('file_and_columns', nargs=-1)
@click.argument('output', nargs=1)
def main(file_and_columns, output):
    '''Merge Multi files by matched column
    '''
    # seperate files and columns
    files = list()
    columns = list()
    for n, each in enumerate(file_and_columns):
        if os.path.isfile(each):
            files.append(each)
        else:
            break
    file_num = len(files)
    if n < file_num - 1:
        columns = file_and_columns[n:]
    # simple check
    if file_num < 2:
        click.secho('File number must more than 2.', color='red')
        return 1
    else:
        if columns:
            if file_num != len(columns):
                click.secho(
                    'File number should equal to column number.', color='red')
                return 1
        else:
            columns = [1 for each in range(file_num)]
        # read files based on merging columns
        df_list = [
            pd.read_table(files[i], index_col=columns[i] - 1)
            for i in range(file_num)
        ]
        # merge files
        merged_df = merge_files(df_list)
        merged_df.to_csv(output, sep='\t')


if __name__ == '__main__':
    main()
