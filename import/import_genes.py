'''
to be used in nextflow
'''
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from models import db_utils
from lib import utils
import argparse
import yaml

def main(params):
    if params.mysql_yaml is None:
        params.mysql_yaml = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mysql.yml')
    with open(params.mysql_yaml, 'rt') as inf:
        mysql_params = yaml.safe_load(inf)
    if params.port:
        mysql_params['mysql']['port'] = params.port
    conn = db_utils.db_connect(mysql_params)
    
    out_header = ['chrom', 'start', 'end', 'ensembl_id', 'symbol', 'mysql_id']

    # import genes
    '''
    no coordinates imported, since the SV files come with gene_ids/symbols
    '''
    cur = conn.cursor()
    genes = utils.get_genes_from_ensembl(params.gtf)
    with open(args.output, 'wt') as outf:
        outf.write('#' + '\t'.join(out_header) + '\n')
        for gene in genes:
            gene['mysql_id'] = db_utils.add_record(cur, 'Gene', gene)
            outf.write('\t'.join([str(gene[h]) for h in out_header]) + '\n')
    conn.commit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', dest='gtf')
    parser.add_argument('--port', dest='port', type=int)
    parser.add_argument('--mysqlyaml', dest='mysql_yaml', default=None)
    parser.add_argument('--output', dest='output', default=None)
    args = parser.parse_args()
    main(args)
