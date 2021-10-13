'''
infer participant id from file name
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import argparse
import yaml


def dictate(row, header, attributes, table):
    # convert data so that they match the datatype in models
    # and 'None', '' to None
    data_type_dict = {
            'int': int,
            'tinyint': lambda x: {'1':1, '0':0, 'True': 1, 'False': 0}[x],
            'float':float
    }
    result = {}
    for item in zip(header, row):
        if item[1] in ('None', ''):
            result[item[0]] = None
        else:
            column = [i for i in attributes[table]['columns'] if i['name'] == item[0]][0]
            column_dt = column['type'].split('(')[0]
            result[item[0]] = data_type_dict.get(column_dt, str)(item[1])

    return result

        
def main(args):
    chunk_size = 1000
    # connect to mysql
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)

    # read table attributes
    FILEPATH = os.path.dirname(os.path.realpath(__file__))
    _default_attributes_file = os.path.join(FILEPATH, '..', 'models', 'attributes.yml')

    attributes = {}

    missing = [i.rstrip() for i in open('missing', 'rt').readlines()]

    with open(_default_attributes_file, 'rt') as inf:
        attributes = yaml.safe_load(inf)

    # write relation tables the last

    #input_files = [os.path.join(args.input_path, f) for f in sorted(os.listdir(args.input_path), key=lambda x: '_' in x)]
    input_files = [os.path.join(args.input_path, f) for f in os.listdir(args.input_path) if f in missing]
    for input_file in input_files:
        print(f"importing {input_file}")
        chrom, tp, table = os.path.basename(input_file).rstrip('.tsv').split('.')
        
        header = []
        records = []
        cur = conn.cursor()
        with open(input_file, 'rt') as inf:
            header = []
            for line in inf:
                row = line.rstrip().split('\t')
                if not header:
                    header = row
                    continue
                row_dict = dictate(row, header, attributes, table)
                records.append(row_dict)
                if len(records) >= chunk_size:
                    db_utils.add_records(cur, table, records)
                    records = []
        if records:
            db_utils.add_records(cur, table, records)

        conn.commit()
        cur.close()

                    
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', dest='input_path', help='the region file generated from the walk')
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    args = parser.parse_args()
    main(args)

