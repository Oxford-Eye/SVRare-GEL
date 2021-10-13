'''
infer participant id from file name
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from lib import Interval_base, utils
import pysam
import argparse
import gzip
import tempfile
import yaml
from collections import Counter
import subprocess

from datetime import datetime

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
                result[item[0]] = item[1].split(',')
        elif item[1] in ('None', ''):
            result[item[0]] = None
        elif item[0] in ('start', 'end', 'CN'):
            result[item[0]] = int(item[1])
        else:
            result[item[0]] = item[1]
    return result

def get_duplicates(intervals, all_intervals, distance_cutoff=0.1):
    # some individuals are sequenced more than once, and those SV files might give the same SV, but with slightly different boundaries
    # if distance is small enough, consider them as the same, and keep a random one that passes the filter
    bad_intervals = set()
    participant_intervals = {}
    for interval in intervals:
        for patient in all_intervals[interval.interval_id]['patients']:
            if patient['patient_id'] not in participant_intervals:
                participant_intervals[ patient['patient_id'] ] = []
            participant_intervals[ patient['patient_id'] ].append({'interval':interval, 'filter': patient['filter'], 'sv_id': patient['sv_id']})

    for p_id, values in participant_intervals.items():
        # sort intervals, PASS in front. So when choose, always choose the ind_i
        these_intervals = sorted(values, key=lambda x: x['filter'] != 'PASS')
        for ind_i in range(len(these_intervals) - 1):
            if (these_intervals[ind_i]['sv_id'], p_id) in bad_intervals:
                continue
            for ind_j in range(ind_i+1, len(these_intervals)):
                if (these_intervals[ind_j]['sv_id'], p_id) in bad_intervals:
                    continue
                this_distance = Interval.get_distance(these_intervals[ind_i]['interval'], these_intervals[ind_j]['interval'])
                if this_distance <= distance_cutoff:
                    bad_intervals.add((these_intervals[ind_j]['sv_id'], p_id))

    return bad_intervals

def tabix_fetch(fname, chrom, start=None, end=None):
    region = chrom
    if start is not None:
        region += f':{start}'
        if end is not None:
            region += f'-{end}'
    cmd = ['tabix', fname, region]
    return subprocess.check_output(cmd).decode('utf8').split('\n')
        
def write_table(out_handles, table, record, attributes):
    if isinstance(record, dict):
        out_handles[table].write('\t'.join([str(record[h['name']]) for h in attributes[table]['columns'] if h['name'] != 'id']) + '\n')
    else:
        out_handles[table].write('\t'.join([str(getattr(record, h['name'])) for h in attributes[table]['columns'] if h['name'] != 'id']) + '\n')

def main(args):

    # get attribultes
    attribute_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'models', 'attributes.yml')
    with open(attribute_file, 'rt') as inf:
        attributes = yaml.safe_load(inf)

    chrom, interval_type, _ = os.path.basename(args.region_file).split('.', 2)
    out_handles = {
            'Interval': open(os.path.join(args.output_path, f"{chrom}.{interval_type}.Interval.tsv"), 'wt'),
            'Interval_Gene': open(os.path.join(args.output_path, f"{chrom}.{interval_type}.Interval_Gene.tsv"), 'wt'),
            'Patient_Interval': open(os.path.join(args.output_path, f"{chrom}.{interval_type}.Patient_Interval.tsv"), 'wt'),
            'Group': open(os.path.join(args.output_path, f"{chrom}.{interval_type}.Group.tsv"), 'wt'),
            'Interval_Group': open(os.path.join(args.output_path, f"{chrom}.{interval_type}.Interval_Group.tsv"), 'wt'),
    }
    # headers for outputs
    for table, outf in out_handles.items():
        outheader = [h['name'] for h in attributes[table]['columns'] if h['name'] != 'id']
        outf.write('\t'.join(outheader) + '\n')

    # read patient table {plate_key: int(patient_id)}
    patient_file_table = {}
    with open(args.patient_file, 'rt') as inf:
        patient_header = []
        for line in inf:
            row = line.rstrip().split('\t')
            if not patient_header:
                patient_header = row
                continue
            row_dict = dict(zip(patient_header, row))
            patient_file_table[row_dict['plate_key']] = int(row_dict['participant_id'])

    with open(args.region_file, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            now = datetime.now().strftime("%H:%M:%S")
            print(now, row[:4])
            start = int(row[1])
            end = int(row[2])
            input_files = row[4].split(',')
            # debug memory leak
            # sizes = sorted([(k,sys.getsizeof(v)) for k,v in list(locals().items())], key=lambda x:x[1], reverse=True) 
            # print(sizes[:10])
            # sizes = sorted([(k,sys.getsizeof(v)) for k,v in list(globals().items())], key=lambda x:x[1], reverse=True) 
            # print(sizes[:10])


            all_intervals = {}
            super_groups = {}
            header = []
            for input_file in input_files:
                # parse participant id and plate key
                if not os.path.isfile(input_file):
                    continue
                plate_key = os.path.basename(input_file).split('.')[0]
                participant_id = patient_file_table[plate_key]
                if not header:
                    with gzip.open(input_file, 'rt') as inf:
                        for line in inf:
                            if line.startswith('#'):
                                header = line.lstrip('#').rstrip().split('\t')
                                break

                tbx = pysam.TabixFile(input_file)
                try:
                    fetcher = tbx.fetch(chrom, start, end)
                except ValueError:
                    print(f"Warn: contig {chrom} not in file {input_file}")
                    continue

                for row in fetcher:
                    if not row:
                        continue
                    row_dict = dictate(row, header)
                    if row_dict['type'] != interval_type:
                        continue
                    row_dict.update({
                        'patient_id': participant_id,
                        'interval_id': '-'.join([chrom, str(row_dict['start']), str(row_dict['end']), interval_type])
                    })
                    # add Interval and Interval_Gene
                    #cur = conn.cursor()
                        # row_dict['interval_id'] = db_utils.add_record(cur, 'Interval', row_dict)

                        # participant_id might be sequenced multiple times. if one of them is pass, pass
                        # 1 = pass, 0 = not pass
                    if row_dict['interval_id'] not in all_intervals:
                        genes = set(row_dict['genes'] + row_dict['genes_at_bnd']) - {''}
                        all_intervals[row_dict['interval_id']] = {
                                'genes': genes,
                                'genes_at_bnd': [i for i in genes if i in row_dict['genes_at_bnd']],
                                'patients': [],
                                'Interval': Interval({
                                    'chrom': chrom,
                                    'start': row_dict['start'],
                                    'end': row_dict['end'],
                                    'type': row_dict['type'],
                                    'interval_id': row_dict['interval_id'],
                                }),
                        }

                    all_intervals[row_dict['interval_id']]['patients'].append({
                        k['name']: row_dict[k['name']] for k in attributes['Patient_Interval']['columns'] if k['name'] in row_dict
                    })
            
            if not all_intervals:
                return
            # write to Interval_Gene
            for interval_id, interval_value in all_intervals.items():
                for gene_id in interval_value['genes']:
                    interval_gene_output = {
                            'interval_id': interval_id,
                            'gene_id': gene_id,
                            'at_bnd': gene_id in interval_value['genes_at_bnd'],
                    }
                    write_table(out_handles, 'Interval_Gene', interval_gene_output, attributes)

            now = datetime.now().strftime("%H:%M:%S")
            print(now, 'before converting to Interval')
            Intervals = [i['Interval'] for i in all_intervals.values()]
            now = datetime.now().strftime("%H:%M:%S")
            print(now, 'before supergrouping')
            intervals_super_group = Interval.group(Intervals)
            super_groups.update(
                {i.mean_id: i for i in intervals_super_group})
            now = datetime.now().strftime("%H:%M:%S")
            print(now, 'after supergrouping')
            # super_group
            for super_group in intervals_super_group:
                #cur = conn.cursor()
                # insert super_group
                super_group.type = 'super'

                # remove duplicates, and insert intervals/interval_gene
                # returns set of (sv_id, patient_id)
                duplicate_intervals = get_duplicates(super_group.intervals, all_intervals)
                now = datetime.now().strftime("%H:%M:%S")
                print(now, 'after a dedup')
                dedup_intervals = []
                for interval in super_group.intervals:
                    for patient in all_intervals[interval.interval_id]['patients']:
                        if (patient['sv_id'], patient['patient_id']) in duplicate_intervals:
                            patient['is_duplicate'] = 1
                        else:
                            patient['is_duplicate'] = 0
                            dedup_intervals.append(patient)
                        write_table(out_handles, 'Patient_Interval', patient, attributes)
                        # write interval
                    write_table(out_handles, 'Interval', interval, attributes)

                        #out_handles['Interval_Gene'].write('\t'.join([str(interval_gene_output[h['name']]) for h in attributes['Interval_Gene']['columns']  if h['name'] != 'id']) + '\n')
                    '''
                    if genes:
                        interval_gene_output = [{
                            'interval_id': interval.interval_id,
                            'gene_id': gene_id,
                            'at_bnd': gene_id in interval.genes_at_bnd,
                            } for gene_id in genes]
                        db_utils.add_records(cur, 'Interval_Gene', interval_gene_output)
                    '''

                # get interval events (Patient_Interval)
                super_group.N_all = len(dedup_intervals)
                super_group.N_pass = len([i for i in dedup_intervals if i['filter'] == 'PASS'])
                for attr in ('stdev_start', 'stdev_end'):
                    stdev = getattr(super_group, attr)
                    if stdev is not None:
                        stdev = float(stdev)
                    setattr(super_group, attr, stdev)

                super_group.group_id = '-'.join(['super', chrom, str(super_group.outer_start), str(super_group.outer_end), interval_type])
                write_table(out_handles, 'Group', super_group, attributes)
                # insert interval_group relation
                for interval in super_group.intervals:
                    outputs = {'group_id': super_group.group_id, 'interval_id': interval.interval_id}
                    write_table(out_handles, 'Interval_Group', outputs, attributes)
                '''
                outputs = [{'group_id': super_group.group_id, 'interval_id': interval.interval_id} for interval in super_group.intervals]
                db_utils.add_records(cur, 'Interval_Group', outputs)
                conn.commit()
                cur.close()
                '''
                # cluster on dedups
                #super_group.intervals = dedup_intervals
                now = datetime.now().strftime("%H:%M:%S")
                print(now, 'before local_cluster')
                local_groups = super_group.cluster()
                now = datetime.now().strftime("%H:%M:%S")
                print(now, 'after local_cluster')
                # add group to intervals
                for group_key, group in local_groups.items():
                    #cur = conn.cursor()
                    group.type = 'local'
                    # get interval events (Patient_Interval)
                    dedup_interval_events = []
                    for interval in group.intervals:
                        for patient in all_intervals[interval.interval_id]['patients']:
                            if (patient['sv_id'], patient['patient_id']) not in duplicate_intervals:
                                dedup_interval_events.append(patient)
                    group.N_all = len(dedup_interval_events)
                    group.N_pass = len([i for i in dedup_interval_events if i['filter'] == 'PASS'])
                    for attr in ('stdev_start', 'stdev_end'):
                        stdev = getattr(group, attr)
                        if stdev is not None:
                            stdev = float(stdev)
                        setattr(group, attr, stdev)
                    group.group_id = '-'.join(['local', chrom, str(group.outer_start), str(group.outer_end), interval_type])
                    write_table(out_handles, 'Group', group, attributes)
                    for interval in group.intervals:
                        outputs = {'group_id':group.group_id, 'interval_id': interval.interval_id}
                        write_table(out_handles, 'Interval_Group', outputs, attributes)
                    '''
                    outputs = [{'group_id': group.group_id, 'interval_id': interval.interval_id} for interval in group.intervals]
                    db_utils.add_records(cur, 'Interval_Group', outputs)
                    conn.commit()
                    cur.close()
                    '''
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--region_file', dest='region_file', help='the region file generated from the walk')
    parser.add_argument('--patient_file', dest='patient_file', help='patient file')
    parser.add_argument('--output_path', dest='output_path', help='output path')
    args = parser.parse_args()
    main(args)
    print('==done==')

