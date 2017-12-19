from collections import Counter
import click


@click.command()
@click.argument('tracking', type=click.File('rb'), required=True)
def main(tracking):
    counts = Counter()
    for eachline in tracking:
        eachline_inf = eachline.strip().split()
        class_code = eachline_inf[3]
        counts.update(class_code)
    for each in counts:
        print '{c}:{n}'.format(c=each, n=counts[each])


if __name__ == '__main__':
    main()
