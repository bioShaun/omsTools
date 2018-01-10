from __future__ import print_function
import click
import subprocess


def check_job(jobname):
    cmd = 'ps ux | grep {j}'.format(j=jobname)
    job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    result = job.communicate()
    return result[0]


@click.command()
@click.argument(
    'jobname',
    type=click.STRING,
    required=True,
)
def main(jobname):
    job_stat = check_job(jobname)
    print(job_stat)


if __name__ == '__main__':
    main()
