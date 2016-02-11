
'''Parse insert_size_metrics file from Picard's CollectInsertSizeMetrics

Prints a single metric to stdout (for easy use with scripts).
Options include:
    MEDIAN_INSERT_SIZE (median)
    MEDIAN_ABSOLUTE_DEVIATION
    MIN_INSERT_SIZE (min)
    MAX_INSERT_SIZE (max)
    MEAN_INSERT_SIZE (mean)
    STANDARD_DEVIATION (sd|stdev|deviation)
    READ_PAIRS
    PAIR_ORIENTATION (orientation)

The metrics file should be in the unfortunate format of:
    ## METRICS CLASS	picard.analysis.InsertSizeMetrics
    <tab-separated metrics names>
    <tab-separated metric values>

Arguments:
    metric file: The insert_size_metrics file.
    metric: A string denoting which metric to return. Default: mean
Return:
    Prints the metric to stdout.
'''

from os.path import isfile
import sys
from textwrap import dedent

usage = dedent('''
    Usage:  python parse_picard_inslen.py FILE [METRIC]

    Synopsis:
        Prints the given metric to standard out.

    Description:
        The metrics file should be in the unfortunate format of:
            ## METRICS CLASS	picard.analysis.InsertSizeMetrics
            <tab-separated metrics names>
            <tab-separated metric values>

    Note:
        The metric and file positions are interchangeable.

    Arguments:
        FILE    The insert_size_metrics file from CollectInsertSizeMetrics.
        METRIC  The string name of the metric. Options include:
                    MEDIAN_INSERT_SIZE (median)
                    MEDIAN_ABSOLUTE_DEVIATION
                    MIN_INSERT_SIZE (min)
                    MAX_INSERT_SIZE (max)
                    MEAN_INSERT_SIZE (mean)
                    STANDARD_DEVIATION (sd|stdev|deviation)
                    READ_PAIRS
                    PAIR_ORIENTATION (orientation)
        -h      Print this message and exit.
''')

if __name__ == '__main__':

    # Metric lookup table
    metric_lookup = {
        'median': 'MEDIAN_INSERT_SIZE',
        'median_dev': 'MEDIAN_ABSOLUTE_DEVIATION',
        'min_ins': 'MIN_INSERT_SIZE',
        'max_ins': 'MAX_INSERT_SIZE',
        'mean_ins': 'MEAN_INSERT_SIZE',
        'stdev': 'STANDARD_DEVIATION',
        'read_pairs': 'READ_PAIRS',
        'orientation': 'PAIR_ORIENTATION',
    }
    metric_synonyms = {
        'sd': 'stdev',
        'deviation': 'stdev',
        'max': 'max_ins',
        'min': 'min_ins',
        'mean': 'mean_ins',
    }

    # help!
    try:
        help_index = sys.argv.index('-h')
        print(usage)
        quit()
    except ValueError:
        pass  # no help request

    # expect two arguments
    try:
        metric_file, metric = sys.argv[1:]
    except ValueError:
        metric_file = sys.argv[1]
        metric = 'mean_ins'

    # allow file to come first!
    if isfile(metric):
        metric_file, metric = [metric, metric_file]

    # lookup actual metric name
    # print usage and raise ValueError if unknown metric
    try:
        # match to a known synonym
        metric = metric_synonyms[metric]
    except KeyError:
        # otherwise, just ensure the metric is lowercase to match the keys
        metric = metric.lower()
    finally:
        try:
            metric = metric_lookup[metric]
        except KeyError:
            print(usage)
            raise ValueError('Unknown metric ' + metric)

    # open the file and find the metric!
    try:
        with open(metric_file, 'r') as fh:
            for line in fh:

                if line.startswith('## Metrics class'.upper()):
                    metric_names = fh.next().split()
                    metric_values = fh.next().split()
                    break

                if metric in line:
                    break

    except IOError:
        print(usage)
        raise

    # parse for desired metric & print (or raise)
    try:
        index = metric_names.index(metric)
        print(metric_values[index])
    except UnboundLocalError:
        print(usage)
        raise ValueError('Desired metric not found: ' + metric)
