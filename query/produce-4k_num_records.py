import sys
import os
sys.path.append('/nas/weka.gel.zone/re_gecip/hearing_and_sight/JingYu/SV_mysql_import')
from models import db_utils
import argparse
import yaml
import pandas as pd

def main(args):
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)
    family_ids = []
    result = []
    for f in os.listdir(args.input_dir):
        report = os.path.join(args.input_dir, f, 'report.csv')
        if not os.path.isfile(report):
            continue
        family_ids.append(f)
        N = 0
        with open(report, 'rt') as inf:
            header = []
            for line in inf:
                if not header:
                    header = line
                    continue
                N += 1
        result.append({
            'family': f,
            'count': N
        })
    columns = (
        'rare_diseases_family_id',
        'normalised_specific_disease',
        'genome_build',
        'participant_ethnic_category',
    )
    place_holder = ','.join(['%s'] * len(family_ids))
    sql = f"SELECT {','.join(columns)} FROM Patient WHERE participant_type = 'Proband' AND rare_diseases_family_id in ({place_holder})"
    cur = conn.cursor()
    cur.execute(sql, family_ids)
    families = {}
    for row in cur.fetchall():
        row_dict = dict(zip(columns, row))
        families[row_dict['rare_diseases_family_id']] = row_dict

    outheader = (
        'family',
        'count',
        'genome_build',
        'normalised_specific_disease',
        'participant_ethnic_category',
    )


    with open(args.output, 'wt') as outf:
        outf.write('\t'.join(outheader) + '\n')
        for row in result:
            if row['family'] not in families:
                continue
            row.update(families[row['family']])
            outf.write('\t'.join([str(row[h]) for h in outheader]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    parser.add_argument('--output', dest='output', default=None)
    parser.add_argument('--input_dir', dest='input_dir', default="report/4k_families_2021-08-08", help='home dir, used to generate igv batch script')
    args = parser.parse_args()
    main(args)
