import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import sqlite3
import argparse
import yaml

def main(args):
    # import patients
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)
    patients =  []
    cur = conn.cursor()
    with open(args.input, 'rt') as inf:
        patients_header = []
        for line in inf:
            row = line.rstrip().split('\t')
            if not patients_header:
                patients_header = row
                continue
            patients.append(dict(zip(patients_header, row)))
    chunk_size = 1000
    for i in range(0, len(patients), chunk_size):
        db_utils.add_records(cur, 'Patient', patients[i:i + chunk_size])
        conn.commit()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='the patient data file')
    parser.add_argument('--port', dest='port', type=int)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    args = parser.parse_args()
    main(args)
