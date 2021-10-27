'''
report per participant (or family_id, rather)
HPO is not reliable as many probands are missing from it. Therefore will use 'normalised_specific_disease' instead.
'''
import sys
import os
sys.path.append('..')
from models import db_utils
from lib import Interval_base, Group
import argparse
import yaml
from collections import defaultdict
import attr
import pysam
import json
from datetime import date
import csv
import gzip

@attr.s
class Family(object):
    proband = attr.ib()
    disease = attr.ib()
    members = attr.ib(factory=list)
    @property
    def member_pids(self):
        if getattr(self, '_member_pids', None) is None:
            self._member_pids = set([i['participant_id'] for i in self.members])
        return self._member_pids

class Interval(Interval_base):
    def __init__(self, row_dict):
        super().__init__(row_dict['chrom'], row_dict['start'], row_dict['end'])
        for k,v in row_dict.items():
            if getattr(self, k, None) is None:
                setattr(self, k, v)

@attr.s
class Report(object):
    conn = attr.ib()
    family_id = attr.ib()
    tbx_gtf = attr.ib()
    outdir = attr.ib()
    distance_cutoff = attr.ib(default=0.5)
    N_else = attr.ib(default=2)
    N_max_genes = attr.ib(default=10)
    igv_padding = attr.ib(default=0.2) # fraction of interval size
    panel_path = attr.ib(default="/nas/weka.gel.zone/re_gecip/hearing_and_sight/JingYu/SV_mysql_import/panelApp")
    out_columns_base = attr.ib(
        default=attr.Factory(lambda:[
            'chrom',
            'start',
            'end',
            'size',
            'type',
            'sv_id',
            'filter',
            'interesting_family_share',
            'interesting_disease_share',
            'interesting_panel_gene',
            'interesting_panel_gene_cross_cds',
            'interesting_HPO_gene',
            'interesting_HPO_gene_cross_cds',
            'N_all',
            'N_share_disease',
            'N_else',
            'genes',
            'genes_cross_cds',
            'other_carriers',
        ])
    )
    out_columns = attr.ib(init=False, default=list)

    interval_columns = attr.ib(
        default=attr.Factory(lambda:[
            'ip.interval_id', # interval_id. maybe there's an elegant way of doing this?
            'ip.chrom',
            'start',
            'end',
            'ip.type',
            'filter',
            'sv_id',
            'genotype',
            'CN',
        ])
    )

    query_columns = attr.ib(
        default=attr.Factory(lambda:[
            'chrom',
            'start',
            'end',
            'interval_id',
            'sv_id',
            'patient_id',
            'normalised_specific_disease',
            'rare_diseases_family_id',
            'genotype',
            'CN'
        ])
    )

    @property
    def panel_genes(self):
        if getattr(self, '_panel_genes', None) is None:
            # not sure if the disease is there
            self._panel_genes = set()
            disease_file = os.path.join(self.panel_path, f"{self.family.disease}.tsv")
            if os.path.isfile(disease_file):
                header = []
                with open(disease_file, 'rt') as inf:
                    for line in inf:
                        row = line.rstrip().split('\t')
                        if not header:
                            header = row
                            continue
                        row_dict = dict(zip(header, row))
                        self._panel_genes.add(row_dict['EnsemblId(GRch38)'])
        return self._panel_genes
    @property
    def beautified_out_columns(self):
        if getattr(self, '_boc', None) is None:
            _boc = []
            self.family
            for column in self.out_columns:
                if column in self.family.member_pids:
                    member = [i for i in self.family.members if i['participant_id'] == column][0]
                    _boc.append(
                        f"{member['relation_to_proband']}(participant_id:{column},disease:{member['disease']})"
                    )
                else:
                    _boc.append(column)
            self._boc = _boc
        return self._boc
        

    @property
    def family(self):
        if getattr(self, '_family', None) is None:
            self._family = self._get_family()
            # put other_carriers at the end. It's so annoying to have it at the front
            self.out_columns = self.out_columns_base + \
                sorted([i['participant_id'] for i in self._family.members], key=lambda x: x!=self._family.proband)
            if 'other_carriers' in self.out_columns_base:
                self.out_columns.remove('other_carriers')
                self.out_columns.append('other_carriers')
        return self._family

    @property
    def HPO_genes(self):
        if getattr(self, '_HPO_genes', None) is None:
            self._HPO_genes = self._get_HPO_genes()
        return self._HPO_genes

    def report(self):
        igv_genome_locations = {
            'GRCh38': '/public_data_resources/IGV/hg38/hg38_local.genome',
            'GRCh37': '/public_data_resources/IGV/hg19/hg19_local.genome',
        }
        os.makedirs(self.outdir, exist_ok=True)
        if not self.family.proband:
            err_file = os.path.join(self.outdir, 'ERROR')
            with open(err_file, 'wt') as outf:
                outf.write("Proband's genome not available!")
            return

        genome_build = [i['genome_build'] for i in self.family.members if i['participant_id'] == self.family.proband][0]
        if genome_build == 'N/A':
            err_file = os.path.join(self.outdir, 'ERROR')
            with open(err_file, 'wt') as outf:
                outf.write("Proband's genome not available!")
            return

        today = date.today().strftime("%Y-%m-%d")
        meta_file = os.path.join(self.outdir, 'meta.yml')
        meta_header = ('distance_cutoff', 'N_else', 'N_max_genes')
        with open(meta_file, 'wt') as outf:
            outf.write(f"time:{today}\n")
            for h in meta_header:
                outf.write(f"{h}:{getattr(self, h)}\n")
        
        report_file = os.path.join(self.outdir, 'report.csv')
        # write report
        igv = []
        with open(report_file, 'wt') as outf:
            csvwriter = csv.writer(outf)
            csvwriter.writerow(self.beautified_out_columns)
            for row in self._get_rare_sv():
                # report
                out_row = [str(row.get(h, '')) for h in self.out_columns]
                csvwriter.writerow(out_row)

                # igv
                size = row['end'] - row['start']
                igv_entry = {
                    'chrom': row['chrom'],
                    'start': row['start'] - int(self.igv_padding * size),
                    'end': row['end'] + int(self.igv_padding * size),
                    'size': size,
                    'sv_id': row['sv_id'],
                }
                igv.append(igv_entry)
        # write igv plot
        igv_file = os.path.join(self.outdir, 'igv_batch.txt')

        igv_plot_dir = os.path.join(self.outdir, 'igv_plots')
        os.makedirs(igv_plot_dir, exist_ok=True)
        with open(igv_file, 'wt') as outf:
            build = [i['genome_build'] for i in self.family.members if i['participant_id'] == self.family.proband][0]
            bam_files_to_load = [i['bam_path'] for i in sorted(self.family.members, key = lambda x: x['participant_id'] != self.family.proband) if i['genome_build'] == build]
            if bam_files_to_load:
                outf.write("new\n")
                outf.write("preference SAM.MAX_VISIBLE_RANGE 1000 \n")
                outf.write(f"genome {igv_genome_locations[build]}\n")
                outf.write(f"snapshotDirectory {igv_plot_dir}\n")
                for bam_file in bam_files_to_load:
                    outf.write(f"load {bam_file}\n")
                for row in igv:
                    chrom, start, end = '', '', ''
                    if build == 'GRCh38':
                        chrom = row['chrom']
                        start = row['start']
                        end = row['end']
                    elif build == 'GRCh37':
                        # need to locate record in one of the vcfs
                        proband_vcf = [i['sv_path'] for i in self.family.members if i['participant_id'] == self.family.proband][0]
                        with gzip.open(proband_vcf, 'rt') as inf:
                            for line in inf:
                                if line.startswith('#'):
                                    continue
                                vcf_row = line.rstrip().split('\t')
                                if vcf_row[2] == row['sv_id']:
                                    chrom, start = vcf_row[0], int(vcf_row[1]) - int(self.igv_padding * row['size'])
                                    for field in vcf_row[7].split(';'):
                                        if field.startswith('END='):
                                            end = int(field.split('=')[1]) + int(self.igv_padding * row['size'])
                                    break

                    outf.write(f"goto {chrom}:{start}-{end}\n")
                    outf.write("sort base\n")
                    outf.write("squish\n")
                    outf.write(f"snapshot {row['chrom']}_{row['start']}-{row['end']}.png\n")
            outf.write('exit')
                    



    def _get_rare_sv(self):

        interval_sql = f"SELECT {','.join(self.interval_columns)} FROM (SELECT group_id FROM `Group` WHERE N_all < 100) G JOIN Interval_Group on G.group_id = Interval_Group.group_id JOIN (SELECT `Interval`.interval_id, chrom, type, start, end, count(*) as c from `Interval` JOIN Patient_Interval ON `Interval`.interval_id = Patient_Interval.interval_id GROUP BY `Interval`.interval_id having c < 50) ip ON ip.interval_id = Interval_Group.interval_id JOIN (SELECT interval_id as ii, filter, genotype, CN, sv_id FROM Patient_Interval WHERE patient_id = %s and is_duplicate = 0) pi ON ip.interval_id = ii"
        seen = set()
        interval_cur = self.conn.cursor()
        interval_cur.execute(interval_sql, (self.family.proband,))
        # should be about 300
        N = 0
        for interval_row in interval_cur.fetchall():
            interval = dict(zip(self.interval_columns, interval_row))
            interval['interval_id'] = interval['ip.interval_id']
            interval['chrom'] = interval['ip.chrom']
            interval['type'] = interval['ip.type']
            if interval['interval_id'] in seen:
                continue

            interval['size'] = interval['end'] - interval['start']
            # find similar intervals
            sql = f"""SELECT {','.join(self.query_columns)} FROM 
                            (SELECT interval_id as ii, chrom, start, end FROM `Interval` WHERE chrom = %s
                                AND type = %s
                                AND NOT (
                                    end < %s
                                    OR end > %s
                                    OR start > %s
                                    OR start < %s
                                )
                                AND (least(end, %s) - greatest(start, %s)) / (greatest(end, %s) - least(start, %s)) > %s
                            ) fil 
                            JOIN Patient_Interval ON Patient_Interval.interval_id = fil.ii
                            JOIN Patient ON patient_id = participant_id
                            WHERE is_duplicate = 0
                    """
            core_cur = self.conn.cursor()
            core_cur.execute(sql, (
                interval['chrom'],                                          # (SELECT * FROM `Interval` WHERE chrom = %s
                interval['type'],                                           #    AND type = %s
                                                                            #    AND NOT
                interval['end'] - interval['size'] * self.distance_cutoff,  #        end < %s
                interval['end'] + 
                    interval['size'] * (1 / (1-self.distance_cutoff) - 1),  #        OR end > %s
                interval['start'] + interval['size']*self.distance_cutoff,  #        OR start > %s
                interval['start'] - 
                    interval['size'] * (1 / (1-self.distance_cutoff) - 1),  #        OR start < %s)
                interval['end'],                                            #   AND (least(end, %s)
                interval['start'],                                          #       - greateset(start, %s) / 
                interval['end'],                                            #       (greatest(end, %s)
                interval['start'],                                          #       - least(start, %s)
                1 - self.distance_cutoff,                                   #   > %s
            ))
            interval_groups = defaultdict(list)
            from collections import Counter
            C = Counter()
            for core_row in core_cur.fetchall():
                query_interval = dict(zip(self.query_columns, core_row))
                C[query_interval['interval_id']] += 1
                if query_interval['patient_id'] == self.family.proband:
                    seen.add(query_interval['interval_id'])
                interval_groups[query_interval['patient_id']].append(query_interval)

            core_cur.close()
            interval['Ns'] = self._get_Ns(interval_groups)
            if interval['Ns']['N_else'] <= self.N_else:
                out_row = self._get_out_row(interval, interval_groups)
                yield out_row
        interval_cur.close()

    def _get_Ns(self, interval_groups):
        # Fail = True, Not fail = False
        N_else = 0
        N_share_disease = 0
        for p_id in interval_groups:
            if p_id not in self.family.member_pids:
                if interval_groups[p_id][0]['normalised_specific_disease'] != self.family.disease:
                    N_else += 1
                else:
                    N_share_disease += 1
        return {'N_else': N_else, 'N_share_disease': N_share_disease}

    def _get_out_row(self, interval, interval_group):
        result = {k:v for k,v in interval.items()}
        result.update({
            'N_all': len(interval_group),
            'N_else': interval['Ns']['N_else'],
            'N_share_disease': interval['Ns']['N_share_disease'],
            'other_carriers': [],
            'genes': [],
        })
        for p_id, value in interval_group.items():
            if value[0]['genotype'] is not None:
                result[p_id] = value[0]['genotype']
            else:
                result[p_id] = f"CN:{value[0]['CN']}"
            if p_id == self.family.proband:
                continue
            if p_id in self.family.member_pids:
                result['interesting_family_share'] = True
            else:
                if value[0]['normalised_specific_disease'] == self.family.disease:
                    result['interesting_disease_share'] = True
                result['other_carriers'].extend(value)
        # dump json for other_carriers
        result['other_carriers'] = json.dumps(result['other_carriers']) if result['other_carriers'] else ''
        
        # for CDS and Genes, treat INV differently (only interested at INV boundaries)
        # cds
        if interval['type'] == 'INV':
            cds = self._get_overlap_cds({
                'chrom': interval['chrom'],
                'start': interval['start'],
                'end': interval['start']+1,
            })
            cds.update(self._get_overlap_cds({
                'chrom': interval['chrom'],
                'start': interval['end'],
                'end': interval['end'] + 1,
            }))
        else:
            cds = self._get_overlap_cds(interval)
        # genes
        columns = ('gene_id', 'symbol')
        if interval['type'] == 'INV':
            sql = f"SELECT {','.join(columns)} FROM Interval_Gene JOIN Gene ON gene_id = ensembl_id WHERE interval_id = %s AND at_bnd = 1"
        else:
            sql = f"SELECT {','.join(columns)} FROM Interval_Gene JOIN Gene ON gene_id = ensembl_id WHERE interval_id = %s"
        cur = self.conn.cursor()
        cur.execute(sql, (interval['interval_id'],))
        for row in cur.fetchall():
            row_dict = dict(zip(columns, row))
            result['genes'].append(row_dict['symbol'])
            if row_dict['gene_id'] in cds:
                result['genes_cross_cds'] = True
                if row_dict['gene_id'] in self.HPO_genes:
                    result['interesting_HPO_gene_cross_cds'] = True
                if row_dict['gene_id'] in self.panel_genes:
                    result['interesting_panel_gene_cross_cds'] = True
            if row_dict['gene_id'] in self.HPO_genes:
                result['interesting_HPO_gene'] = True
            if row_dict['gene_id'] in self.panel_genes:
                result['interesting_panel_gene'] = True

        if len(result['genes']) > self.N_max_genes:
            result['genes'] = result['genes'][:self.N_max_genes] + ['...']
        result['genes'] = ','.join(result['genes'])
        cur.close()
        return result
            

    def _get_family(self):
        result = {
            'members': [],
            'proband': None,
            'disease': None
        }
        # family member
        columns = (
            'participant_id',
            'normalised_specific_disease',
            'biological_relationship_to_proband',
            'participant_type',
            'plate_key',
            'path',
            'genome_build',
        )
        sql = f"SELECT {','.join(columns)} FROM Patient WHERE rare_diseases_family_id = %s"
        cur= self.conn.cursor()
        cur.execute(sql, (self.family_id,))
        for row in cur.fetchall():
            row_dict = dict(zip(columns, row))
            bam_path = os.path.join(row_dict['path'], row_dict['plate_key'], 'Assembly', f"{row_dict['plate_key']}.bam")
            sv_path = os.path.join(row_dict['path'], row_dict['plate_key'], 'Variations', f"{row_dict['plate_key']}.SV.vcf.gz")
            if row_dict['participant_type'] == 'Proband':
                result['disease'] = row_dict['normalised_specific_disease']
                result['proband'] = row_dict['participant_id']
                result['members'].append({
                    'participant_id': row_dict['participant_id'],
                    'relation_to_proband': 'proband',
                    'disease': row_dict['normalised_specific_disease'],
                    'genome_build': row_dict['genome_build'],
                    'bam_path': bam_path,
                    'sv_path': sv_path,
                })
            else:
                result['members'].append({
                    'participant_id': row_dict['participant_id'],
                    'relation_to_proband': row_dict['biological_relationship_to_proband'],
                    'disease': row_dict['normalised_specific_disease'],
                    'genome_build': row_dict['genome_build'],
                    'bam_path': bam_path,
                    'sv_path': sv_path,
                })
        return Family(**result)
    
    def _get_overlap_cds(self, interval):
        '''
        what cds an interval crosses
        tbx_gtf is a pysam.TabixFile
        '''
        # note that the gtf contigs don't have chr.
        result = {}
        chrom = interval['chrom'].lstrip('chr')
        if chrom == 'M':
            chrom = 'MT'
        for line in self.tbx_gtf.fetch(chrom, interval['start'], interval['end']):
            row = line.rstrip().split('\t')
            if row[2] != 'CDS':
                continue
            info = {}
            for field in row[-1].split('; '):
                key, value = field.split(' ',1)
                info[key] = value.split('"')[1]
            if info['gene_id'] not in result:
                result[info['gene_id']] = {}
            if info['transcript_id'] not in result[info['gene_id']]:
                result[info['gene_id']][info['transcript_id']] = []
            result[info['gene_id']][info['transcript_id']].append(int(info['exon_number']))
        return result

    def _get_HPO_genes(self):
        # get interesting genes according to proband's HPO.
        # if HPO is missing, return empty set
        sql = "SELECT gene_id FROM (SELECT * FROM Patient_HPO WHERE patient_id = %s) ph JOIN HPO_Gene ON ph.hpo_id = HPO_Gene.hpo_id"
        genes = set()
        cur = self.conn.cursor()
        cur.execute(sql, (self.family.proband,))
        for row in cur.fetchall():
            genes.add(row[0])
        cur.close()
        return genes

def main(args):
    # connect to mysql
    with open(args.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if args.port:
        mysql_params['mysql']['port'] = args.port
    conn = db_utils.db_connect(mysql_params)
    if args.outdir is None:
        args.outdir = os.path.join('report', args.family_id)
    
    report = Report(conn, args.family_id, pysam.TabixFile('/re_gecip/hearing_and_sight/JingYu/SV_mysql_query/data/Homo_sapiens.GRCh38.98.sorted.gtf.gz'), args.outdir, args.interval['distance_cutoff'], args.interval['N_else'])
    report.report()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--family_id', dest='family_id', help='the family id is usually the same as the proband id')
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
    print('==done==')


