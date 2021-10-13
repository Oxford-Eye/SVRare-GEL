import pysam
import argparse
import os
import bisect
import gzip

def get_header(fname):
    with gzip.open(fname, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                return line.lstrip('#').rstrip().split('\t')
def merge_intervals(intervals):
    '''
    input: intervals = {(start, end): set([files]), ...}
    output: result same structure as input 
    '''

    merged = []
    for interval in sorted(list(intervals.keys()), key=lambda x: x[0]):
        if not merged or merged[-1][1] < interval[0]:
            merged.append([interval[0], interval[1], intervals[interval]])
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
            merged[-1] = [merged[-1][0], max(merged[-1][1], interval[1]), merged[-1][2] | intervals[interval]  ]
    return {(i[0],i[1]):i[2] for i in merged}
def main(args):
    intervals = {}
    header = []

    for fname in os.listdir(args.input_path):
        if not fname.endswith('.tsv.gz'):
            continue
        input_file = os.path.join(args.input_path, fname)
        if not header:
            header = get_header(input_file)
        
        tbx = pysam.TabixFile(input_file)
        try:
            fetcher = tbx.fetch(args.chrom)
        except ValueError:
            print(f"Warn: contig {args.chrom} not in file {input_file}")
            continue
        for row in tbx.fetch(args.chrom):
            row_dict = dict(zip(header, row.split('\t')))
            if row_dict['type'] != args.type:
                continue
            this = (int(row_dict['start']), int(row_dict['end']))
            if this not in intervals:
                intervals[this] = set()
            intervals[this].add(input_file)
    merged = merge_intervals(intervals)

            
    with open(args.output, 'wt') as outf:
        header = ['#chrom', 'start', 'end', 'count', 'files']
        outf.write('\t'.join(header) + '\n')
        for interval in merged:
            out_row = [str(i) for i in [
                args.chrom,
                interval[0],
                interval[1],
                len(merged[interval]),
                ','.join(merged[interval])
            ]]
            outf.write('\t'.join(out_row) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', dest='chrom', help='which chrom to process')
    parser.add_argument('--input_path', dest='input_path', help='input files separated by space')
    parser.add_argument('--type', dest='type' )
    parser.add_argument('--output', dest='output' )

    args = parser.parse_args()
    main(args)

