import gzip
import os
import sys
import subprocess
import tempfile
import gzip



def liftover_file(input_bed, chain_file):
    '''
    simple wrapper of liftOver
    make sure to add a row id at the fourth column
    '''
    tmp_output_fd, tmp_output = tempfile.mkstemp()
    tmp_unmapped_fd, tmp_unmapped = tempfile.mkstemp()
    subprocess.run(['liftOver', input_bed, chain_file, tmp_output, tmp_unmapped], check=True)
    with open(tmp_output, 'rt') as inf:
        for line in inf:
            yield line
    os.close(tmp_output_fd)
    os.close(tmp_unmapped_fd)
    os.remove(tmp_output)
    os.remove(tmp_unmapped)

def makebed(vcf, bed):
    with gzip.open(vcf, 'rt') as inf, open(bed, 'wt') as outf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            chrom = row[0]
            start = stop = row[1]
            row_id = row[2]
            for field in row[7].split(';'):
                if field.startswith('END='):
                    stop = field.split('=')[1]
            if int(start) > int(stop):
                # might be MantaINV. will check in later code.
                start, stop = stop, start
            out_row = [chrom, start, stop, row_id]
            outf.write('\t'.join(out_row) + '\n')

def get_genes_from_ensembl(ensembl_file):
    result = []
    with gzip.open(ensembl_file, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if row[2] != 'gene': continue
            chrom = row[0]
            # change MT to M
            if chrom == 'MT':
                chrom = 'M'
            this_record = {
                    'chrom': f"chr{chrom}",
                    'start': int(row[3]),
                    'end': int(row[4]),
                    }
            for field in row[-1].split('; '):
                key, val = field.split(' ')
                this_record[key] = val.split('"')[1]
            this_record['ensembl_id'] = this_record.pop('gene_id')
            this_record['symbol'] = this_record.pop('gene_name')
            result.append(this_record)
    return result
