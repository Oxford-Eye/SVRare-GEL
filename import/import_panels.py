import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import argparse
import yaml
import csv
import re

def main(args):
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)
    tables = ('Panel','Patient_Panel')
    sqls = db_utils.create_tables(tables=tables)
    cur = conn.cursor()
    for sql in sqls:
        cur.execute(sql)
    conn.commit()
    cur.close()
    # panels
    for f in os.listdir(args.panel_path):
        panel_name = os.path.splitext(f)[0]
        print(f)
        outputs = []
        with open(os.path.join(args.panel_path, f), 'rt') as inf:
            header = []
            for line in inf:
                row = line.rstrip().split('\t')
                if not header:
                    header = row
                    continue
                row_dict = dict(zip(header, row))
                moi = ''
                if row_dict['Model_Of_Inheritance']:
                    moi = re.split('[^A-Za-z-]', row_dict['Model_Of_Inheritance'])[0]

                try:
                    gel_status = int(row_dict['GEL_Status'])
                except ValueError:
                    gel_status = 0

                output = {
                    'gene_id': row_dict['EnsemblId(GRch38)'],
                    'symbol': row_dict['Gene Symbol'],
                    'name': panel_name,
                    'moi': moi,
                    'status': gel_status,
                }
                outputs.append(output)
        cur = conn.cursor()
        db_utils.add_records(cur, 'Panel', outputs)
        conn.commit()
        cur.close()
    # patient_panel
    with open(args.patient_panel, 'rt') as inf:
        csvreader = csv.reader(inf)
        header = []
        cur = conn.cursor()
        for line in csvreader:
            if not header:
                header = line
                continue
            row_dict = dict(zip(header, line))
            output = {
                'patient_id': int(row_dict['Participant Id']),
                'panel_name': row_dict['Panel Name']
            }
            # sometimes panel name doesn't exist
            if os.path.isfile(os.path.join(args.panel_path, f"{output['panel_name']}.tsv")):
                sql = 'SELECT participant_id FROM Patient WHERE participant_id = %s'
                cur.execute(sql, (output['patient_id'],))
                if list(cur.fetchall()):
                    try:
                        db_utils.add_record(cur, 'Patient_Panel', output)
                    except:
                        print(output)
                        raise
        conn.commit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    parser.add_argument('--panel_path', dest='panel_path', default="panelApp")
    parser.add_argument('--patient_panel', dest='patient_panel', default="input/panels_applied_2021-08-09_13-15-59.csv")
    args = parser.parse_args()
    main(args)


