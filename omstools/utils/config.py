#!/usr/bin/env python

""" omsTools config """

import pkg_resources
import click

# click mutually exclusive option class


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


# config for click help
CLICK_CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

# gtf tools
# GTF tools includes class, functions and dicts from assemblyline
gtf_tools = dict()
for entry_point in pkg_resources.iter_entry_points('gtf'):
    gtf_tools.update({entry_point.name: entry_point.load()})
