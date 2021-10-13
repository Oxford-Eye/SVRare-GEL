'''
infer participant id from file name
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
from lib import Interval_base, utils
import pysam
import argparse
import gzip
import tempfile
import yaml
from collections import Counter
import subprocess

def construct_dict_key(row_dict):
    return ':'.join([row_dict['chrom'], str(row_dict['start']), str(row_dict['end']), row_dict['type']])

class Interval(Interval_base):
    def __init__(self, row_dict):
        for k,v in row_dict.items():
            setattr(self, k, v)
        self.size = self.end - self.start

def dictate(row, header):
    # try to convert str to number if possible
    # and 'None' to None
    result = {}
    for item in zip(header, row.split('\t')):
        if item[0] in ('genes', 'genes_at_bnd'):
            if item[1] == '':
                result[item[0]] = []
            else:
                result[item[0]] = [int(i) for i in item[1].split(',')]
        elif item[1] in ('None', ''):
            result[item[0]] = None
        elif item[0] in ('start', 'end', 'CN'):
            result[item[0]] = int(item[1])
        else:
            result[item[0]] = item[1]
    return result

def get_duplicates(intervals, distance_cutoff=0.1):
    # some individuals are sequenced more than once, and those SV files might give the same SV, but with slightly different boundaries
    # if distance is small enough, consider them as the same, and keep a random one that passes the filter
    bad_intervals = set()
    participant_intervals = {}
    for interval in intervals:
        if interval.patient_id not in participant_intervals:
            participant_intervals[ interval.patient_id ] = []
        participant_intervals[ interval.patient_id ].append(interval)

    for p_id, values in participant_intervals.items():
        # sort intervals, PASS in front. So when choose, always choose the ind_i
        these_intervals = sorted(values, key=lambda x: x.filter != 'PASS')
        for ind_i in range(len(these_intervals) - 1):
            if (these_intervals[ind_i].sv_id, these_intervals[ind_i].patient_id) in bad_intervals:
                continue
            for ind_j in range(ind_i+1, len(these_intervals)):
                if (these_intervals[ind_j].sv_id, these_intervals[ind_j].patient_id) in bad_intervals:
                    continue
                this_distance = Interval.get_distance(these_intervals[ind_i], these_intervals[ind_j])
                if this_distance <= distance_cutoff:
                    bad_intervals.add((these_intervals[ind_j].sv_id, these_intervals[ind_j].patient_id))

    return bad_intervals

def tabix_fetch(fname, chrom, start=None, end=None):
    region = chrom
    if start is not None:
        region += f':{start}'
        if end is not None:
            region += f'-{end}'
    cmd = ['tabix', fname, region]
    return subprocess.check_output(cmd).decode('utf8').split('\n')
        
def main(args):
    # connect to mysql
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)

    # get columns of `Interval`
    # patient_interval_columns = db_utils.get_columns(conn, 'Patient_Interval')

    all_intervals = {'LOSS':[], 'INV':[], 'GAIN':[]}
    interval_gene_output = []
    super_groups = {'LOSS':{}, 'INV':{}, 'GAIN':{}}
    header = []
    input_files = (os.path.join(args.input_path, f) for f in os.listdir(args.input_path) if f.endswith('.tsv.gz'))
    for input_file in input_files:
        # parse participant id and plate key
        participant_id, plate_key = os.path.basename(input_file).split('.')[0].split('_', 1)
        participant_id = int(participant_id)
        if not header:
            with gzip.open(input_file, 'rt') as inf:
                for line in inf:
                    if line.startswith('#'):
                        header = line.lstrip('#').rstrip().split('\t')
                        break

        # sometimes it gets OSError for getting index for unknown reason.
        try:
            tbx = pysam.TabixFile(input_file)
            try:
                fetcher = tbx.fetch(args.chrom)
            except ValueError:
                print(f"Warn: contig {args.chrom} not in file {input_file}")
                continue
        except OSError:
            print(f"pysam failed. try subprocess")
            fetcher = tabix_fetch(input_file, args.chrom)

        for row in fetcher:
            if not row:
                continue
            row_dict = dictate(row, header)
            # key = construct_dict_key(row_dict)
            # add Interval_Patient 
            row_dict.update({
                'patient_id': participant_id,
            })
            # add Interval and Interval_Gene
            #cur = conn.cursor()
                # row_dict['interval_id'] = db_utils.add_record(cur, 'Interval', row_dict)

                # participant_id might be sequenced multiple times. if one of them is pass, pass
                # 1 = pass, 0 = not pass
            if row_dict['type'] not in all_intervals:
                # only deals with the three types
                continue
            all_intervals[row_dict['type']].append(row_dict)
            # insert patient_interval
            #db_utils.add_record(cur, 'Patient_Interval', all_intervals[row_dict['type']][key])

            #conn.commit()
            #cur.close()
            
    for tp in ('LOSS', 'GAIN', 'INV'):
        intervals = [Interval(i) for i in all_intervals[tp]]
        if not intervals:
            continue
        intervals_super_group = Interval.group(intervals)
        super_groups[tp].update(
            {i.mean_id: i for i in intervals_super_group})
        # super_group
        for super_group in intervals_super_group:
            cur = conn.cursor()
            # insert super_group
            super_group.type = 'super'

            # remove duplicates, and insert intervals/interval_gene
            # returns set of (sv_id, patient_id)
            duplicate_intervals = get_duplicates(super_group.intervals)
            dedup_intervals = []
            for interval in super_group.intervals:
                if (interval.sv_id, interval.patient_id) in duplicate_intervals:
                    interval.is_duplicate = 1
                else:
                    interval.is_duplicate = 0
                    dedup_intervals.append(interval)
                # insert interval
                interval.interval_id = db_utils.add_record(cur, 'Interval', interval)
                # insert patient_interval
                db_utils.add_record(cur, 'Patient_Interval', interval)
                # insert interval_gene
                genes = set(interval.genes + interval.genes_at_bnd) - {''}
                if genes:
                    interval_gene_output.extend([{
                        'interval_id': interval.interval_id,
                        'gene_id': gene_id,
                        'at_bnd': gene_id in interval.genes_at_bnd,
                        } for gene_id in genes])
                    db_utils.add_records(cur, 'Interval_Gene', interval_gene_output)

            # get interval events (Patient_Interval)
            super_group.N_all = len(dedup_intervals)
            super_group.N_pass = len([i for i in dedup_intervals if i.filter == 'PASS'])
            for attr in ('stdev_start', 'stdev_end'):
                stdev = getattr(super_group, attr)
                if stdev is not None:
                    stdev = float(stdev)
                setattr(super_group, attr, stdev)
            try:
                super_group_id = db_utils.add_record(cur, 'Group', super_group)
            except:
                print(args.chrom, tp, [i.sv_id for i in super_group.intervals])
                raise
            # insert interval_group relation
            outputs = [{'group_id': super_group_id, 'interval_id': interval.interval_id} for interval in super_group.intervals]
            db_utils.add_records(cur, 'Interval_Group', outputs)
            conn.commit()
            cur.close()
            # cluster on dedups
            super_group.intervals = dedup_intervals
            local_groups = super_group.cluster()
            # add group to intervals
            for group_key, group in local_groups.items():
                cur = conn.cursor()
                group.type = 'local'
                # get interval events (Patient_Interval)
                group.N_all = len(group.intervals)
                group.N_pass = len([i for i in group.intervals if i.filter == 'PASS'])
                for attr in ('stdev_start', 'stdev_end'):
                    stdev = getattr(super_group, attr)
                    if stdev is not None:
                        stdev = float(stdev)
                    setattr(group, attr, stdev)
                group_id = db_utils.add_record(cur, 'Group', group)
                outputs = [{'group_id': group_id, 'interval_id': interval.interval_id} for interval in group.intervals]
                db_utils.add_records(cur, 'Interval_Group', outputs)
                conn.commit()
                cur.close()
                    
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', dest='chrom', help='which chrom to process')
    parser.add_argument('--input_path', dest='input_path', help='input files separated by space')
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    args = parser.parse_args()
    main(args)

