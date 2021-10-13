import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import argparse
import yaml
from lib import Interval_base, utils
from collections import Counter

def get_phenotypes(input_file):
    return

def main(args):
    # connect to mysql
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)
    
    # read phenotype file
    phenotypes = get_phenotypes(args.phenotype_file)

    sql = "SELECT `Interval`.interval_id, `Interval`.chrom, `Interval`.start, `Interval`.end FROM `Interval` JOIN Interval_Gene ON `Interval`.interval_id = Interval_Gene.interval_id JOIN Gene ON Interval_Gene.gene_id = Gene.ensembl_id WHERE Gene.symbol = %s AND `Interval`.type = 'LOSS'"
    cur = conn.cursor()
    cur.execute(sql, ('OPN4',))
    intervals = []
    for row in cur.fetchall():
        interval = Interval_base(row[1],row[2],row[3])
        interval.interval_id = row[0]
        intervals.append(interval)
        sql = "SELECT interval_id,patient_id,filter,source,genotype,CN,rare_diseases_family_id,normalised_specific_disease FROM Patient_Interval JOIN Patient ON patient_id = participant_id where interval_id = %s and is_duplicate = 0"
        curr = conn.cursor()
        curr.execute(sql, (interval.interval_id,))
        interval.patients = []
        interval.neighbours = []
        for patient_row in curr.fetchall():
            patient_row_dict = dict(zip(['interval_id','patient_id', 'filter', 'source', 'genotype', 'CN', 'f_id', 'disease'], patient_row))
            interval.patients.append(patient_row_dict)
        curr.close()
    intervals.sort()
    intervals_count = Counter({i.interval_id:len(i.patients) for i in intervals})
    print(intervals[0])

    for i in range(len(intervals) - 1):
        for j in range(i+1, len(intervals)):
            if Interval_base.get_distance(intervals[i], intervals[j]) < args.interval['distance_cutoff']:
                intervals[i].neighbours.extend(intervals[j].patients)
                intervals[j].neighbours.extend(intervals[i].patients)
            else:
                break

    for interval in intervals:
        # same patient count as 1
        #count = len(interval.patients) + len(interval.neighbours)
        count = len(set([i['patient_id'] for i in interval.patients + interval.neighbours]))
        if  count > args.interval['N_all']:
            continue
        print(interval.interval_id, count)
        print('patients')
        print(interval.patients)
        #print([i for i in interval.patients if i['genotype'] == 'HOM' or i['CN'] == 0])
        print('neighbours')
        print(interval.neighbours)
        #print(len(interval.neighbours))
        print('====')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--config', dest='config', default='query_config.yml')
    args = parser.parse_args()
    with open(args.config, 'rt') as inf:
        config = yaml.safe_load(inf)
        for key, val in config.items():
            setattr(args, key, val)
    main(args)

