'''
parse dbvar
return an Interval instance, and a dictionary for later annotation 
'''
import pysam
from lib import Interval_base


'''
clinical_assertion can be one or more of these values, delimited by ';'
===
Benign
Benign/Likely benign
Conflicting interpretations of pathogenicity
Likely benign
Likely pathogenic
Pathogenic
Pathogenic/Likely pathogenic
Uncertain significance
association
clinical_assertion
conflicting data from submitters
drug response
not provided
protective
risk factor
===
'''
class Dbvar(object):
    tbx_header = [
        'chrom',
        'start',
        'end',
        'sample_count',
        'sub_type',
        'method',
        'analysis',
        'platform',
        'study',
        'variant',
        'clinical_assertion',
        'clinvar_assertion',
        'bin_size'
    ]
    annotation_header = [
        'id',
        'cnv_type',
        'sub_type',
        'distance',
        'sample_count',
        'method',
        'analysis',
        'platform',
        'study',
        'clinical_assertion',
        'clinvar_assertion',
        'bin_size'
    ]

    def __init__(self, fname, cnv_type):
        self.fname = fname
        self.cnv_type = cnv_type
        self.tbx = pysam.TabixFile(fname)

    def get_dbvar(self, chrom, start, end):
        result = []
        try:
            it = self.tbx.fetch(chrom, max(0, start-1), end)
        except ValueError:
            return result
        for line in it:
            row_dict = dict(zip(self.tbx_header, line.split('\t')))
            # clinical_assertion sometimes can be too long. limit it to no more than 20 items.
            clinical_assertion = row_dict.get(
                'clinical_assertion', '').split(';')
            if len(clinical_assertion) > 20:
                row_dict['clinical_assertion'] = ';'.join(
                    clinical_assertion[:21]) + ' ...'
            interval = Interval_base(
                row_dict['chrom'],
                int(row_dict['start']),
                int(row_dict['end'])
            )
            interval.annotation_header = self.annotation_header
            interval.cnv_type = self.cnv_type
            for key, val in row_dict.items():
                if getattr(interval, key, None) is None:
                    setattr(interval, key, val)
            interval.distance = interval.get_distance(
                interval, Interval_base(chrom, start, end))
            result.append(interval)
        return result
