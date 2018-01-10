from __future__ import print_function
import click
import subprocess
import os

PS_HEADER = 'USER       \
PID \
%CPU \
%MEM    \
VSZ   \
RSS \
TTY      \
STAT \
START   \
TIME \
COMMAND'

INFO_BAR = '=' * 72
EMPTY_LINE = "\n"

CONFIRM_STRING = {
    'kill': 'Are you sure you want to kill these jobs?',
    'rename': 'Wrong job name?',
}


def abort_if_false(ctx, param, value):
    if not value:
        ctx.abort()


class LINUXJOB(object):

    def __init__(self, job_name):
        self.job_name = job_name
        self.job_stat = None
        self.job_pids = list()

    def _sys_run(self, cmd):
        job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        result = job.communicate()
        return result

    def _sys_clear(self):
        os.system('clear')

    def _sys_kill(self, pid):
        cmd = 'kill -9 {pid}'.format(pid=pid)
        self._sys_run(cmd)

    def get_job_id(self):
        if self.job_stat is not None:
            job_lines = self.job_stat.split('\n')
            for each in job_lines:
                if each.strip() != '':
                    self.job_pids.append(each.split()[1])

    def check_job(self):
        cmd = 'ps ux | grep {t.job_name} | grep -v "grep {t.job_name}\|oms_kill\|killjob\|ps ux"'.format(t=self)
        self.job_stat = self._sys_run(cmd)[0]
        return self.job_stat

    def show_job(self):
        self.check_job()
        click.secho(EMPTY_LINE)
        click.secho('Jobs named with [{t.job_name}]:'.format(t=self),
                    fg='blue')
        click.secho(INFO_BAR, fg='blue')
        click.secho(PS_HEADER)
        click.secho(self.job_stat)
        click.secho(INFO_BAR, fg='blue')

    def _confirm_job(self):
        self._sys_clear()
        self.show_job()
        if click.confirm(CONFIRM_STRING['kill']):
            return True
        else:
            click.secho(INFO_BAR, fg='blue')
            if click.confirm(CONFIRM_STRING['rename']):
                click.secho(INFO_BAR, fg='blue')
                self.job_name = click.prompt('Please enter a name', type=str)
                click.secho(INFO_BAR, fg='blue')
                click.secho('New job to kill is [{t.job_name}]'.format(t=self))
                return self._confirm_job()
            else:
                return False

    def kill_job(self):
        if self._confirm_job():
            click.secho(INFO_BAR, fg='blue')
            self.get_job_id()
            if self.job_pids:
                for each_pid in self.job_pids:
                    click.secho('Kill job [{p}]'.format(p=each_pid),
                                fg='red')
                    self._sys_kill(each_pid)
            else:
                click.secho('Nothing to kill!', fg='green')
        else:
            click.secho(INFO_BAR, fg='blue')
            click.secho('Not kill anything!', fg='green')


@click.command()
@click.argument(
    'job_name',
    type=click.STRING,
    required=True,
)
def main(job_name):
    '''Kill job by job name in linux
    '''
    my_jobs = LINUXJOB(job_name)
    my_jobs.kill_job()


if __name__ == '__main__':
    main()
