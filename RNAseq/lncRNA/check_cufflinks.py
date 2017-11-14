import click
import pandas as pd
import os


def read_cufflinks_log(log_file, log=True):
    sample_list = []
    log_inf = open(log_file).readlines()
    for n, eachline in enumerate(log_inf):
        if log:
            if 'Assembling transcripts and estimating abundances' in eachline:
                if (log_inf[n + 1]).strip()[-4:] == '100%':
                    if log_inf[n + 2].strip().split()[-1] == 'finished':
                        sample_list.append(log_inf[n + 2].strip().split()[0])
        else:
            if 'finished' in eachline:
                sample = eachline.split()[1].strip('"')
                sample_list.append(sample)
    return sample_list


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
    type=click.Path(dir_okay=False),
    required=True,
    help='cufflinks run status.')
def main(script_dir, output):
    run_stat_list = []
    for root, dirs, files in os.walk(script_dir):
        script_name = os.path.basename(os.path.abspath(root))
        if 'script' in script_name:
            run_name = '{s}_run'.format(s=script_name)
            done_name = '{s}_done'.format(s=script_name)
            run_list = []
            done_list = []
            for each_file in files:
                each_file_path = os.path.join(root, each_file)
                if 'log' in each_file:
                    done_list.extend(read_cufflinks_log(each_file_path))
                elif each_file[-3:] == '.sh':
                    run_list.extend(read_cufflinks_log(each_file_path, False))
            run_list = list(set(run_list))
            done_list = list(set(done_list))
            run_list_value = [-1 for each in run_list]
            done_list_value = [1 for each in done_list]
            run_stat_list.append(
                pd.Series(run_list_value, index=run_list, name=run_name))
            run_stat_list.append(
                pd.Series(done_list_value, index=done_list, name=done_name))
    run_stat_df = pd.concat(run_stat_list, axis=1)
    run_stat_df.loc[:, 'sum'] = run_stat_df.T.sum().T
    run_stat_df.to_csv(output, sep='\t', na_rep=0)


if __name__ == '__main__':
    main()
