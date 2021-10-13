'''
HPO was imported after everything else was done. Use it with care if starting from scratch
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import argparse
import yaml

def main(args):
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)

    cur = conn.cursor()
    tables = ('HPO','HPO_Gene','Patient_HPO')
    sqls = db_utils.create_tables(tables=tables)
    for sql in sqls:
        cur.execute(sql)
    conn.commit()

    patient_HPO_outputs = []
    HPO_outputs = {}
    HPO_Gene_outputs = []
    with open(args.patient_hpo, 'rt') as inf:
        header = []
        for line in inf:
            row = line.rstrip().split('\t')
            if not header:
                header = row
                continue
            row_dict = dict(zip(header, row))
            if row_dict['Hpo Present'] != 'Yes':
                continue
            patient_id = int(row_dict['Participant Id'])

            sql = "SELECT * FROM Patient WHERE participant_id = %s"
            cur.execute(sql,(patient_id,))
            if cur.fetchone() is None:
                continue
            
            output = {
                'hpo_id': int(row_dict['Hpo Id'][3:]),
                'name': row_dict['Hpo Term'],
                'patient_id': patient_id,
            }
            patient_HPO_outputs.append(output)
            HPO_outputs[output['hpo_id']] = output

    
    with open(args.hpo_gene, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            gene = row[3]
            sql = 'SELECT ensembl_id FROM Gene WHERE symbol = %s'
            cur.execute(sql, (gene,))
            query_result = cur.fetchone()
            if query_result is None:
                continue
            gene_id = query_result[0]
            hpo_id = int(row[0][3:])
            if hpo_id not in HPO_outputs:
                HPO_outputs[hpo_id] = {
                    'hpo_id': hpo_id,
                    'name': row[1]
                }
            HPO_Gene_outputs.append({
                'gene_id': gene_id,
                'hpo_id': hpo_id
            })

    for output in HPO_outputs.values():
        db_utils.add_record(cur, 'HPO', output)

    conn.commit()
    for output in HPO_Gene_outputs:
        db_utils.add_record(cur, 'HPO_Gene', output)
    conn.commit()
    for output in patient_HPO_outputs:
        try:
            db_utils.add_record(cur, 'Patient_HPO', output)
        except:
            print(output)
            raise
    conn.commit()

    cur.close()
    conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    parser.add_argument('--hpo_gene', dest='hpo_gene', default=None)
    parser.add_argument('--patient_hpo', dest='patient_hpo', default=None)
    args = parser.parse_args()
    main(args)

