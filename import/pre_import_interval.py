import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import sqlite3
import argparse
from models import db_utils
    # recreate tables
def main(args):
    conn = sqlite3.connect(args.db)
    db_utils.create_table(conn, 'Interval')
    db_utils.create_table(conn, 'Patient_Interval')
    db_utils.create_table(conn, 'Interval_Gene')
    conn.commit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', dest='db', help='database')
    args = parser.parse_args()
    main(args)
