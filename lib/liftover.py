'''
liftover binary is required.
`module load bio/liftOver/1.0
liftover files live in /public_data_resources/liftover/
'''
import subprocess
import tempfile
import os

def liftover_variants(variants, build_from=None, build_to=None, chainfile=None):
    '''
    if chainfile is provided, build_from and build_to will be ignored
    '''
    if chainfile is None:
        if not (build_from is None and build_to is None) or (build_from == 'GRCh37' and build_to == 'GRCh38'):
            errmsg = f"build_from has to be GRCh37 ({build_from} given) and build_to has to be GRCh38 ({build_to} given), or you can provide a chainfile."
            raise ValueError(errmsg)
        chainfile = '/public_data_resources/liftover/grch37_to_grch38.over.chain.gz'
    
    tmp_bed_fd, tmp_bed = tempfile.mkstemp()
    makebed(variants, tmp_bed)
    result = {}
    for line in liftover_file(tmp_bed, chainfile):
        row = line.rstrip().split('\t')
        orig_variant = row[3]
        _,_,ref,alt = orig_variant.split('-')
        result[row[3]] = '-'.join([row[0], row[1], ref, alt])

    os.close(tmp_bed_fd)
    os.remove(tmp_bed)
    return result

def liftover_intervals(intervals, build_from=None, build_to=None, chainfile=None):
    '''
    if chainfile is provided, build_from and build_to will be ignored
    '''
    if chainfile is None:
        if not (build_from is None and build_to is None) or (build_from == 'GRCh37' and build_to == 'GRCh38'):
            errmsg = f"build_from has to be GRCh37 ({build_from} given) and build_to has to be GRCh38 ({build_to} given), or you can provide a chainfile."
            raise ValueError(errmsg)
        chainfile = '/public_data_resources/liftover/grch37_to_grch38.over.chain.gz'
    
    tmp_bed_fd, tmp_bed = tempfile.mkstemp()
    makebed_intervals(intervals, tmp_bed)
    result = {}
    for line in liftover_file(tmp_bed, chainfile):
        row = line.rstrip().split('\t')
        result[row[3]] = {
            'chrom': row[0],
            'start': int(row[1]),
            'end': int(row[2])
        }

    os.close(tmp_bed_fd)
    os.remove(tmp_bed)
    return result

def liftover_vcf(in_vcf, out_vcf, chainfile, ref_index):
    import gzip
    all_header = []
    chrom_header = []
    variant_array = []
    variant_dict = {}
    ref_file = os.path.splitext(ref_index)[0]

    # contigs
    contigs = []
    with open(ref_index, 'rt') as inf:
        for line in inf:
            row = line.rstrip().split('\t')
            contigs.append({
                'chrom': row[0],
                'length': row[1],
            })

    # variants from in_vcf
    with gzip.open(in_vcf, 'rt') as inf:
        for line in inf:
            if line.startswith('##'):
                if not (line.startswith('##reference=') or line.startswith('##contig=')):
                    all_header.append(line)
                continue
            row = line.rstrip().split('\t')
            if line.startswith('#'):
                chrom_header = row
                continue
            row_dict = dict(zip(chrom_header, row))
            variant_id = '-'.join([row_dict['#CHROM'], row_dict['POS'], row_dict['REF'], row_dict['ALT']])
            variant_array.append(variant_id)
            variant_dict[variant_id] = row_dict
    
    tmp_bed_fd, tmp_bed = tempfile.mkstemp()
    makebed(variant_array, tmp_bed)

    # lifted over variants, sort
    result = []
    for line in liftover_file(tmp_bed, chainfile):
        row = line.rstrip().split('\t')
        orig_variant = row[3]
        out_row = variant_dict[orig_variant]
        out_row.update({
            '#CHROM': row[0],
            'POS': row[1],
        })
        result.append(out_row)
    result.sort(key = lambda x: (x['#CHROM'], int(x['POS'])))
    tmp_out_fd, tmp_out = tempfile.mkstemp()
    with open(tmp_out, 'wt') as outf:
        for line in all_header:
            outf.write(line)
        reference_line = f"##reference=file://{ref_file}\n"
        outf.write(reference_line)
        for contig in contigs:
            outline = f"##contig=<ID={contig['chrom']},length={contig['length']}>\n"
            outf.write(outline)
        outf.write('\t'.join(chrom_header) + '\n')
        for out_row in result:
            outf.write('\t'.join([out_row[h] for h in chrom_header]) + '\n')
    
    with open(out_vcf, 'w') as outf:
        subprocess.call(['bgzip', '-c', tmp_out], stdout=outf)
    subprocess.run(['tabix', '-p', 'vcf', out_vcf])

    os.close(tmp_bed_fd)
    os.remove(tmp_bed)
    os.close(tmp_out_fd)
    os.remove(tmp_out)


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

def makebed_intervals(intervals, bed):
    with open(bed, 'wt') as outf:
        for interval in intervals:
            if isinstance(interval, dict):
                out_row = [interval['chrom'], interval['start'], interval['end'], interval['interval_id']]
            else:
                out_row = [interval.chrom, interval.start, interval.end, interval.interval_id]
            outf.write('\t'.join([str(i) for i in out_row]) + '\n')

def makebed(variants, bed):
    with open(bed, 'wt') as outf:
        for variant in variants:
            chrom, pos, ref, alt = variant.split('-')
            out_row = [chrom, pos, str(int(pos)+1), variant]
            outf.write('\t'.join(out_row) + '\n')
