'''
to be used in nextflow
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
import argparse
import yaml

def main(params):
    # make tables
    if params.mysql_yaml is None:
        params.mysql_yaml = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mysql.yml')
    with open(params.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if params.port:
        mysql_params['mysql']['port'] = params.port
    conn = db_utils.db_connect(mysql_params)
    cur = conn.cursor()
    sqls = db_utils.create_tables(tables=params.tables)
    for sql in sqls:
        cur.execute(sql)
    conn.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', dest='port', type=int)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    parser.add_argument('--tables', dest='tables', nargs='+', default=None)
    args = parser.parse_args()
    main(args)
