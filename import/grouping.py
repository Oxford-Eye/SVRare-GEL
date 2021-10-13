import argparse
from lib import Interval_base, utils
from models import db_utils
import sqlite3

CHROMOSOMES = (
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX',
    'chrY',
    'chrM'
)
class Interval(Interval_base):
    def __init__(self, row_dict):
        for k,v in row_dict.items():
            setattr(self, k, v)
        self.size = self.end - self.start

def get_intervals_from_db(conn, chrom=None, start=None, end=None):
    sql = '''SELECT * FROM Interval'''
    if chrom is not None:
        sql += ''' WHERE chrom = ?'''
        if start is not None:
            sql += ''' AND start >= ?'''
        if end is not None:
            sql += ''' AND end <= ?'''
    cur = conn.cursor()
    cur.execute(sql, [i for i in (chrom, start, end) if i is not None])
    for row in cur.fetchall():
        d = {}
        for idx, col in enumerate(cur.description):
            d[col[0]] = row[idx]
        yield Interval(d)

def main(args):
    conn = sqlite3.connect(args.db)
    # remove Group and Interval_Group
    cur = conn.cursor()
    for table_name in ('Group', 'Interval_Group'):
        sql = f"DROP TABLE IF EXISTS `{table_name}`"
        db_utils.create_table(conn, table_name)
    conn.commit()

    for chrom in CHROMOSOMES:
        print(f"doing {chrom}")
        super_groups = {'LOSS': {}, 'GAIN': {}, 'INV': {}}  # core_id: group
        all_intervals = list(get_intervals_from_db(conn, chrom))

        for tp in ('LOSS', 'GAIN', 'INV'):
            intervals = [i for i in all_intervals if i.type ==
                         tp]
            intervals_super_group = Interval.group(intervals)
            super_groups[tp].update(
                {i.mean_id: i for i in intervals_super_group})
            # super_group
            for super_group in intervals_super_group:
                # insert super_group
                super_group.type = 'super'
                # get interval events (Patient_Interval)
                sql = f'''SELECT filter from Patient_Interval WHERE interval_id in ({','.join([str(i.id) for i in super_group.intervals])})'''
                cur.execute(sql)
                N_all = 0
                N_pass = 0
                for interval_event in cur.fetchall():
                    N_all += 1
                    if interval_event[0] == 'PASS':
                        N_pass += 1
                super_group.N_all = N_all
                super_group.N_pass = N_pass
                super_group_id = db_utils.add_record(cur, 'Group', super_group)
                # insert interval_group relation
                outputs = [{'group_id': super_group_id, 'interval_id': interval.id} for interval in super_group.intervals]
                db_utils.add_records(cur, 'Interval_Group', outputs)
                # cluster
                local_groups = super_group.cluster()
                # add group to intervals
                for group_key, group in local_groups.items():
                    group.type = 'local'
                    # get interval events (Patient_Interval)
                    sql = f'''SELECT filter from Patient_Interval WHERE interval_id in ({','.join([str(i.id) for i in super_group.intervals])})'''
                    cur.execute(sql)
                    N_all = 0
                    N_pass = 0
                    for interval_event in cur.fetchall():
                        N_all += 1
                        if interval_event[0] == 'PASS':
                            N_pass += 1
                    group.N_all = N_all
                    group.N_pass = N_pass
                    group_id = db_utils.add_record(cur, 'Group', group)
                    outputs = [{'group_id': group_id, 'interval_id': interval.id} for interval in group.intervals]
                    db_utils.add_records(cur, 'Interval_Group', outputs)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('--chrom', dest='chrom', help='which chrom to group')
    parser.add_argument('--db', dest='db', help='sqlite3 db path')
    args = parser.parse_args()
    main(args)
