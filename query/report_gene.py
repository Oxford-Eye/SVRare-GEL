'''
gene or region, not both
'''
import sys
import os
sys.path.append('..')
from models import db_utils
from lib import Interval_base, Group, liftover_intervals
import argparse
import yaml
from collections import Counter, defaultdict
import attr
import pysam
import json
from datetime import date
import csv

class Interval(Interval_base):
    def __init__(self, row_dict):
        super().__init__(row_dict['chrom'], row_dict['start'], row_dict['end'])
        for k,v in row_dict.items():
            if getattr(self, k, None) is None:
                setattr(self, k, v)

@attr.s
class Report(object):
    conn = attr.ib()
    tbx_gtf = attr.ib()
    outdir = attr.ib()
    gene = attr.ib(default=None)
    region = attr.ib(default=None)
    interval_type = attr.ib(default=None)
    distance_cutoff = attr.ib(default=0.5)
    N_else = attr.ib(default=2)
    N_max_genes = attr.ib(default=0)
    igv_padding = attr.ib(default=0.2)
    track_max = attr.ib(default=5)
    cluster_interval_max = attr.ib(default=50)
    out_columns = attr.ib(
        default=attr.Factory(lambda:[
            'chrom',
            'outer_start',
            'outer_end',
            'outer_size',
            'type',
            'N_all_P',
            'dominant_disease_P',
            'N_disease_P',
            'disease_ratio_P',
            'dominant_HPO_P',
            'N_HPO_P',
            'HPO_ratio_P',
            'N_all_F',
            'dominant_disease_F',
            'N_disease_F',
            'disease_ratio_F',
            'dominant_HPO_F',
            'N_HPO_F',
            'HPO_ratio_F',
            'genes',
            'genes_cross_cds',
            'participants',
        ])
    )
    
    def report(self):
        igv_genome_locations = {
            'GRCh38': '/public_data_resources/IGV/hg38/hg38_local.genome',
            'GRCh37': '/public_data_resources/IGV/hg19/hg19_local.genome',
        }
        os.makedirs(self.outdir, exist_ok=True)
        # meta
        today = date.today().strftime("%Y-%m-%d")
        meta_file = os.path.join(self.outdir, 'meta.yml')
        meta_header = ('distance_cutoff', 'N_else', 'region', 'gene')
        with open(meta_file, 'wt') as outf:
            outf.write(f"time:{today}\n")
            for h in meta_header:
                outf.write(f"{h}:{getattr(self, h)}\n")
        # cluster
        igv_batch = {'GRCh38':{}, 'GRCh37':{}}
        
        igv_outfile = os.path.join(self.outdir, 'igv_batch.txt')
        igv_plot_dir = os.path.join(self.outdir, 'igv_plots')
        report_outfile = os.path.join(self.outdir, 'report.csv')
        with open(report_outfile, 'wt') as outf:
            csvwriter = csv.writer(outf)
            csvwriter.writerow(self.out_columns)
            for tp in ('LOSS', 'GAIN', 'INV'):
                intervals = [i for i in self.intervals if i.type == tp]
                if not intervals:
                    continue
                group = Group(intervals)
                clusters = group.cluster(eps=self.distance_cutoff)
                for cluster in clusters.values():
                    if len(cluster.intervals) >= self.cluster_interval_max:
                        continue
                    out_row = self._get_out_row(cluster)
                    csvwriter.writerow(out_row)
                    # igv
                    # now cluster has participants
                    # use outer coordinates (outer_id)
                    cluster.outer_id = '-'.join([cluster.chrom, str(cluster.outer_start), str(cluster.outer_end), cluster.type])
                    for build in ('GRCh38', 'GRCh37'):
                        igv_intervals = sorted([i['interval'] for i in cluster.participants.values() if i['interval'].genome_build == build], key = lambda x: x.filter != 'PASS')[:self.track_max]
                        if igv_intervals:
                            if build == 'GRCh37':
                                # liftover
                                lifted = liftover_intervals([{
                                    'chrom': cluster.chrom,
                                    'start': cluster.outer_start,
                                    'end': cluster.outer_end,
                                    'interval_id': cluster.outer_id
                                }], chainfile="/public_data_resources/liftover/hg38ToHg19.over.chain.gz")
                                if not lifted:
                                    continue
                                lifted = lifted[cluster.outer_id]
                                cluster_id = '-'.join([lifted['chrom'], str(lifted['start']), str(lifted['end'])])
                                size = lifted['end'] - lifted['start']
                                igv_batch[build][cluster.outer_id] = {
                                    'chrom': lifted['chrom'],
                                    'start': max(lifted['start'] - int(size * self.igv_padding), 1),
                                    'end': lifted['end'] + int(size * self.igv_padding),
                                    'bam_paths': [],
                                }
                            elif build == 'GRCh38':
                                size = cluster.outer_end - cluster.outer_start
                                igv_batch[build][cluster.outer_id] = {
                                    'chrom': cluster.chrom,
                                    'start': max(cluster.outer_start - int(size * self.igv_padding), 1),
                                    'end': cluster.outer_end + int(size * self.igv_padding),
                                    'bam_paths': [],
                                 }
                            for interval in igv_intervals:
                                bam_path = os.path.join(interval.path, interval.plate_key, 'Assembly', f"{interval.plate_key}.bam")
                                igv_batch[build][cluster.outer_id]['bam_paths'].append(bam_path)
            with open(igv_outfile, 'wt') as outf:
                for build in ('GRCh38', 'GRCh37'):
                    if not igv_batch[build]:
                        continue
                    for interval, value in igv_batch[build].items():
                        outf.write('new\n')
                        outf.write('preference SAM.MAX_VISIBLE_RANGE 1000\n')
                        outf.write(f"genome {igv_genome_locations[build]}\n")
                        outf.write(f"snapshotDirectory {igv_plot_dir}\n")
                        for bam_path in value['bam_paths']:
                            outf.write(f"load {bam_path}\n")
                        coordinate = f"{value['chrom']}:{value['start']}-{value['end']}"
                        outf.write(f"goto {coordinate}\n")
                        outf.write("sort base\n")
                        outf.write("squish\n")
                        outf.write(f"snapshot {interval}-{build}.png\n")
                outf.write('exit\n')

    def _get_out_row(self, cluster):
        # type
        cluster.type = cluster.intervals[0].type
        cluster.mean_size = cluster.mean_end - cluster.mean_start
        cluster.outer_size = cluster.outer_end - cluster.outer_start
        # genes
        genes = set()
        for interval in cluster.intervals:
            genes.update(interval.symbols)
        cluster.genes = ','.join(genes)
        # for CDS and Genes, treat INV differently (only interested at INV boundaries)
        # cds
        if cluster.type == 'INV':
            cds = self._get_overlap_cds({
                'chrom': cluster.chrom,
                'start': cluster.mean_start,
                'end': cluster.mean_start+1,
            })
            cds.update(self._get_overlap_cds({
                'chrom': cluster.chrom,
                'start': cluster.mean_end,
                'end': cluster.mean_end+1,
            }))
        else:
            cds = self._get_overlap_cds({
                'chrom': cluster.chrom,
                'start': cluster.mean_start,
                'end': cluster.mean_end,
            })
        cluster.genes_cross_cds = ','.join(sorted(list(cds.keys()), key=lambda x: x != self.gene) )
        # participants
        participants = {}
        families = {}
        # put PASS ones in the last so they will come and replace those that do not pass
        for interval in sorted(cluster.intervals, key = lambda x: x.filter == 'PASS'):
            if interval.participant_id not in participants:
                participants[interval.participant_id] = {
                    'interval': interval,
                    'HPO': set(interval.hpos),
                    'disease': interval.normalised_specific_disease,
                    'genotype': interval.genotype,
                    'CN': interval.CN,
                    'family_id': interval.family_id,
                }
            else:
                participants[interval.participant_id]['interval'] = interval
            if interval.family_id not in families:
                families[interval.family_id] = {
                    'interval': interval,
                    'HPO': set(interval.hpos)
                }
            else:
                families[interval.family_id]['interval'] = interval
                families[interval.family_id]['HPO'].update(interval.hpos)
        cluster.participants = participants
        cluster.N_all_P = len(participants)
        disease_P_counter = Counter()
        HPO_P_counter = Counter()
        for interval in participants.values():
            disease_P_counter[interval['interval'].normalised_specific_disease] += 1
            HPO_P_counter.update(interval['HPO'])
        cluster.dominant_disease_P, cluster.N_disease_P = disease_P_counter.most_common(1)[0]
        cluster.dominant_HPO_P, cluster.N_HPO_P = HPO_P_counter.most_common(1)[0]
        cluster.disease_ratio_P = cluster.N_disease_P / cluster.N_all_P
        cluster.HPO_ratio_P = cluster.N_HPO_P / cluster.N_all_P

        cluster.N_all_F = len(families)
        disease_F_counter = Counter()
        HPO_F_counter = Counter()
        for interval in families.values():
            disease_F_counter[interval['interval'].normalised_specific_disease] += 1
            HPO_F_counter.update(interval['HPO'])
        cluster.dominant_disease_F, cluster.N_disease_F = disease_F_counter.most_common(1)[0]
        cluster.dominant_HPO_F, cluster.N_HPO_F = HPO_F_counter.most_common(1)[0]
        cluster.disease_ratio_F = cluster.N_disease_F / cluster.N_all_F
        cluster.HPO_ratio_F = cluster.N_HPO_F / cluster.N_all_F

        return [getattr(cluster, h) for h in self.out_columns]

    @property
    def intervals(self):
        if getattr(self, '_intervals', None) is None:
            self._intervals = list(self._get_intervals())
        return self._intervals
    def _get_intervals(self):
        columns = (
            '`Interval`.chrom as chrom',
            '`Interval`.start as start',
            '`Interval`.end as end',
            '`Interval`.type as type',
            '`Interval`.interval_id as interval_id',
            'participant_id',
            'sv_id',
            'filter',
            'genotype',
            'CN',
            'rare_diseases_family_id as family_id',
            'biological_relationship_to_proband',
            'normalised_specific_disease',
            'participant_type',
            'plate_key',
            'path',
            'genome_build',
        )
        '''
        sql = f"""SELECT {','.join(columns)} FROM `Interval` 
            JOIN Interval_Group ON `Interval`.interval_id = Interval_Group.interval_id
            JOIN `Group` ON Interval_Group.group_id = `Group`.group_id
            JOIN Patient_Interval ON Patient_Interval.interval_id = `Interval`.interval_id
            JOIN Patient ON patient_id = participant_id
            LEFT JOIN Patient_HPO ON participant_id = Patient_HPO.patient_id
            LEFT JOIN HPO on Patient_HPO.hpo_id = HPO.hpo_id
            WHERE N_all < 100  AND is_duplicate = 0
            """
        '''
        params = []
        sql = f"""SELECT {','.join(columns)} FROM (
            SELECT * FROM `Group` WHERE N_all < 100"""
        if self.region:
            chrom, rest = self.region.split(':')
            start, end = [int(i) for i in rest.split('-')]
            sql += " AND CHROM = %s AND outer_start <= %s AND outer_end >= %s"
            params.extend([chrom, end, start])
        sql += """) G 
            JOIN Interval_Group ON G.group_id = Interval_Group.group_id
            JOIN `Interval` ON `Interval`.interval_id = Interval_Group.interval_id
            JOIN Patient_Interval ON `Interval`.interval_id = Patient_Interval.interval_id
            JOIN Patient ON patient_id = participant_id
            """
        if self.gene:
            sql += """ JOIN Interval_Gene ON `Interval`.interval_id = Interval_Gene.interval_id
                JOIN (SELECT * FROM Gene WHERE symbol = %s) gene ON gene_id = ensembl_id
                """
            params.append(self.gene)
        if args.interval_type:
            sql += """ AND `Interval`.type = %s """
            params.append(self.interval_type)
            
        sql += ' ORDER BY `Interval`.start, `Interval`.end'
        cur = self.conn.cursor()
        cur.execute(sql, params)
        processed_columns = process_columns(columns)
        # result 
        result = {}
        patients = set()
        interval_ids = set()
        interval_count = Counter()
        row_dicts = []
        for row in cur.fetchall():
            row_dict = dict(zip(processed_columns, row))
            row_dicts.append(row_dict)
            interval_count[row_dict['interval_id']] += 1
        # remove intervals with too popular counts
        for row_dict in row_dicts:
            if interval_count[row_dict['interval_id']] >=100:
                continue
            unique_id = '-'.join([row_dict['interval_id'], str(row_dict['participant_id'])])
            if unique_id not in result:
                result[unique_id] = Interval(row_dict)
                patients.add(row_dict['participant_id'])
                interval_ids.add(row_dict['interval_id'])
        if not result:
            return []
        #HPO
        hpo_sql = f"""SELECT patient_id, GROUP_CONCAT(DISTINCT name SEPARATOR ',') FROM Patient_HPO
            JOIN HPO ON Patient_HPO.hpo_id = HPO.hpo_id
            WHERE patient_id IN ({','.join(['%s'] * len(patients))}) GROUP BY patient_id"""
        hpo_cur = self.conn.cursor()
        hpo_cur.execute(hpo_sql, list(patients))
        hpos = {}
        for row in hpo_cur.fetchall():
            hpos[row[0]] = row[1].split(',')
        hpo_cur.close()
        #genes
        genes = {}
        if self.region:
            where_in = ','.join(['%s'] * len(interval_ids))
            gene_sql = f"""SELECT interval_id, GROUP_CONCAT(DISTINCT symbol SEPARATOR ',') FROM Interval_Gene
                JOIN Gene ON Interval_Gene.gene_id = ensembl_id
                WHERE interval_id IN ({where_in}) GROUP BY interval_id"""
            gene_cur = self.conn.cursor()
            gene_cur.execute(gene_sql, list(interval_ids))
            for row in gene_cur.fetchall():
                genes[row[0]] = row[1].split(',')
            gene_cur.close()

        for v in result.values():
            v.hpos = hpos.get(v.participant_id, ['N/A'])
            if self.gene:
                v.symbols = [self.gene]
            else:
                v.symbols = genes.get(v.interval_id, [])
        return result.values()

    def _get_overlap_cds(self, interval):
        '''
        what cds an interval crosses
        tbx_gtf is a pysam.TabixFile
        '''
        # note that the gtf contigs don't have chr.
        result = {}
        for line in self.tbx_gtf.fetch(interval['chrom'].lstrip('chr'), interval['start'], interval['end']):
            row = line.rstrip().split('\t')
            if row[2] != 'CDS':
                continue
            info = {}
            for field in row[-1].split('; '):
                key, value = field.split(' ',1)
                info[key] = value.split('"')[1]
            if info['gene_name'] not in result:
                result[info['gene_name']] = {}
            if info['transcript_id'] not in result[info['gene_name']]:
                result[info['gene_name']][info['transcript_id']] = []
            result[info['gene_name']][info['transcript_id']].append(int(info['exon_number']))
        return result

def process_columns(raw_columns):
    return [i.split('as ')[-1] for i in raw_columns]

def main(args):
    if getattr(args, 'gene', None) and getattr(args, 'region', None):
        msg = 'Only provide one of --gene and --region, not both.'
        raise ValueError(msg)
    # connect to mysql
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)

    report = Report(
        conn,
        pysam.TabixFile('/re_gecip/hearing_and_sight/JingYu/SV_mysql_query/data/Homo_sapiens.GRCh38.98.sorted.gtf.gz'),
        args.outdir,
        args.gene,
        args.region,
        args.interval_type,
        args.interval['distance_cutoff'],
        args.interval['N_else'],
    )
    report.report()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene', dest='gene', help='gene symbol')
    parser.add_argument('--region', dest='region', help='chrom:start-end')
    parser.add_argument('--interval_type', dest='interval_type', help='[LOSS|GAIN|INV]', default=None)
    parser.add_argument('--port', dest='port', type=int, default=None)
    parser.add_argument('--config', dest='config', default='query_config.yml')
    parser.add_argument('--outdir', dest='outdir', default=None)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    args = parser.parse_args()
    with open(args.config, 'rt') as inf:
        config = yaml.safe_load(inf)
        for key, val in config.items():
            setattr(args, key, val)
    main(args)



