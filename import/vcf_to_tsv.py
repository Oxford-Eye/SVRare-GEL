'''
convert GRCh37 to GRCh38 if necessary
types:
    INV, LOSS, GAIN
change into tsv
columns:
    chrom
    start
    end
    sv_id
    filter
    source
    genotype
    CN
    type
    genes
'''

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from lib import Interval_base, utils
import gzip
import argparse
import tempfile
import pysam
import subprocess

def get_genes_from_line(row_dict, build, genes, csqt_header=None):
    '''
    GRCh37 gives ensembl_gene_id
    GRCh38 gives symbol under CSQT
    will return
        {value: set(), type: ensembl_id|symbol}
    '''
    result = set()
    for info in row_dict['INFO'].split(';'):
        if build == 'GRCh37' and info.startswith('ensembl_gene_id'):
            ensembl_ids = set( info.split('=')[1].split(',') )
            result.update([i['ensembl_id'] for i in genes if i['ensembl_id'] in ensembl_ids])
        elif build == 'GRCh38' and info.startswith('CSQT'):
            symbols = set()
            if csqt_header is None:
                raise ValueError('Need to provide csqt_header to parse for GRCh38')
            for csqt_field in info.split('=')[1].split(','):
                csqt = dict(zip(csqt_header, csqt_field.split('|')))
                symbols.add(csqt['HGNC'])
            result.update([i['ensembl_id'] for i in genes if i['symbol'] in symbols])
    return result

class Interval(Interval_base):
    def __init__(self, row_dict):
        ID = row_dict['ID'].split(':')
        chrom = row_dict.get('chrom', row_dict['CHROM'])
        start = int(row_dict['start'])
        self.info = {}
        for field in row_dict['INFO'].split(';'):
            if '=' in field:
                key,val = field.split('=')
                self.info[key] = val
            else:
                self.info[field] = True

        # if liftover, row_dict will have an end entry. if not, get it from info
        end = int(row_dict.get('end', self.info['END']))

        super().__init__(chrom, start, end)
        self.type = None
        if ID[0] == 'Canvas':
            self.source = 'Canvas'
            self.type = ID[1]
        else:
            # only deal with MantaDEL and MantaINV
            self.source = 'Manta'
            if ID[0] == 'MantaDEL':
                self.type = 'LOSS'
            elif ID[0] == 'MantaINV':
                self.type = 'INV'
        self.alt = row_dict['ALT']
        self.sv_id = row_dict['ID']
        self.filter = row_dict['FILTER']
        self.super_groups = []
        self.groups = []
        genotype_dict = dict(zip(row_dict['FORMAT'].split(
            ':'), row_dict['genotype'].split(':')))
        if 'GT' not in genotype_dict:
            self.genotype = None
        else:
            genotype = genotype_dict['GT'].split('/')
            if '.' in genotype:
                self.genotype = genotype
            elif len(set(genotype)) == 1:
                self.genotype = 'HOM'
            else:
                self.genotype = 'HET'
        if 'CN' not in genotype_dict:
            self.CN = None
        else:
            self.CN = int(genotype_dict['CN'])

def get_genes_at_bnd(chrom, pos, gtf_tbx, gtf_header, args):

    bnd_start = max(1, pos - args.padding)
    bnd_end = pos + args.padding
    result = set()
    for row in gtf_tbx.fetch(chrom, bnd_start, bnd_end):
        row_dict = dict(zip(gtf_header, row.split('\t')))
        result.add(row_dict['ensembl_id'])
    return result


def main(args):
    '''
    Files in GRCh37 is annotated differently to those in GRCh38
        all_annotated_intervals_dict = {}
    Also GRCh37 needs to be lifted over
    only deals with Canvas, and MantaDEL/MantaINV
    '''

    vcf_header = ['CHROM', 'start', 'ID', 'REF', 'ALT',
                  'QUAL', 'FILTER', 'INFO', 'FORMAT', 'genotype']
    out_header = ['chrom', 'start', 'end', 'sv_id', 'filter', 'source', 'genotype', 'CN', 'type', 'genes', 'genes_at_bnd']

    # read contigs and gtf header
    contigs = set(subprocess.check_output(['tabix', '-l', args.gtf]).decode('utf8').split('\n'))
    gtf_tbx = pysam.TabixFile(args.gtf)
    gtf_header = []
    with gzip.open(args.gtf, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                gtf_header = line.lstrip('#').rstrip().split('\t')
                break
    # read gtf_genes
    gtf_genes = []
    with gzip.open(args.gtf, 'rt') as inf:
        header = []
        for line in inf:
            row = line.lstrip('#').rstrip().split()
            if not header:
                header = row
                continue
            gtf_genes.append(dict(zip(header, row)))
    

    # GRCh37? get new coordinates
    if args.build == 'GRCh37':
        liftover_dict = {}
        fd, bed37 = tempfile.mkstemp()
        utils.makebed(args.input, bed37)
        for line in utils.liftover_file(bed37, args.liftOver_chainfile):
            row = line.rstrip().split('\t')
            liftover_dict[row[3]] = {
                'chrom': row[0],
                'start': row[1],
                'end': row[2],
            }
        os.close(fd)
        os.remove(bed37)
    # get CSQT header
    csqt_header = None
    if args.build == 'GRCh38':
        with gzip.open(args.input, 'rt') as inf:
            for line in inf:
                if line.startswith('##INFO=<ID=CSQT'):
                    csqt_header = line.split('"')[1].split('Format: ')[1].split('|')

    with gzip.open(args.input, 'rt') as inf:
        print('#' + '\t'.join(out_header))
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            row_dict = dict(zip(vcf_header, row))
            ID = row_dict['ID'].split(':')
            if row_dict['ALT'] == '.':
                continue
            if not (ID[0].startswith('Canvas') or ID[0].startswith('MantaDEL') or ID[0].startswith('MantaINV')):
                continue
            if args.pass_only and row_dict['FILTER'] != 'PASS':
                continue
            if args.build == 'GRCh37':
                if row_dict['ID'] not in liftover_dict:
                    # liftover failed for this record
                    continue
                row_dict.update(liftover_dict[row_dict['ID']])
            genes = ','.join([str(i) for i in get_genes_from_line(row_dict, args.build, gtf_genes, csqt_header)])
            interval = Interval(row_dict)
            interval.genes = genes
            genes_at_bnd = set()
            if interval.chrom in contigs:
                for pos in (interval.start, interval.end):
                    genes_at_bnd.update(get_genes_at_bnd(interval.chrom, pos, gtf_tbx, gtf_header, args))
            interval.genes_at_bnd = ','.join([str(i) for i in genes_at_bnd])
            # insert_interval will also add table rowid to interval.
            if interval.type is None:
                continue
            print('\t'.join([str(getattr(interval, h)) for h in out_header]))

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='path of input file')
    parser.add_argument('--build', dest='build', help='build [GRCh37 | GRCh38]')
    parser.add_argument('--liftOver_chainfile', dest='liftOver_chainfile', help='path to the liftover chain file')
    parser.add_argument('--pass_only', dest='pass_only', default=False, action='store_true',  help='get pass only intervals?')
    parser.add_argument('--padding', dest='padding', type=int, default=100)
    parser.add_argument('--gtf', dest='gtf', help='human ensembl gtf in GRCh38 (sorted and tabixed and with mysql id)')
    args = parser.parse_args()
    main(args)

