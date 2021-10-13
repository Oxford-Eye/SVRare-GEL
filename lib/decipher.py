'''
parse decipher
return an Interval instance, and a dictionary for later annotation
'''
import pysam
from lib import Interval_base


class Decipher(object):
    tbx_header = [
        'population_cnv_id',
        'chrom',
        'start',
        'end',
        'deletion_observations',
        'deletion_frequency',
        'deletion_standard_error',
        'duplication_observations',
        'duplication_frequency',
        'duplication_standard_error',
        'observations',
        'frequency',
        'standard_error',
        'type',
        'sample_size',
        'study',
    ]
    annotation_header = [
        'id',
        'cnv_type',
        'distance',
        'sample_count',
        'frequency',
        'standard_error',
        'type',
        'sample_size',
        'study',
    ]

    def __init__(self, fname, cnv_type):
        self.fname = fname
        self.tbx = pysam.TabixFile(fname)
        self.cnv_type = cnv_type

    def get_decipher(self, chrom, start, end):
        result = []
        try:
            it = self.tbx.fetch(chrom, max(0, start-1), end)
        except ValueError:
            return result

        for line in it:
            row_dict = dict(zip(self.tbx_header, line.split('\t')))
            for key in ('deletion_observations', 'duplication_observations', 'observations', 'sample_size'):
                row_dict[key] = int(row_dict[key])
            for key in ('deletion_standard_error', 'duplication_standard_error', 'deletion_frequency', 'duplication_frequency', 'frequency', 'standard_error'):
                row_dict[key] = float(row_dict[key])
            # note that each row might have both deletions and duplications

            if row_dict['deletion_observations'] > 0 and self.cnv_type == 'LOSS':
                interval = Interval_base(
                    row_dict['chrom'],
                    int(row_dict['start']),
                    int(row_dict['end'])
                )
                interval.annotation_header = self.annotation_header
                interval.cnv_type = 'LOSS'
                interval.sample_size = row_dict['sample_size']
                interval.sample_count = row_dict['deletion_observations']
                interval.frequency = row_dict['deletion_frequency']
                interval.standard_error = row_dict['deletion_standard_error']
                interval.type = row_dict['type']
                interval.study = row_dict['study']
                interval.distance = interval.get_distance(
                    interval, Interval_base(chrom, start, end))
                result.append(interval)
            if row_dict['duplication_observations'] > 0 and self.cnv_type == 'GAIN':
                interval = Interval_base(
                    row_dict['chrom'],
                    int(row_dict['start']),
                    int(row_dict['end'])
                )
                interval.annotation_header = self.annotation_header
                interval.cnv_type = 'GAIN'
                interval.sample_size = row_dict['sample_size']
                interval.sample_count = row_dict['duplication_observations']
                interval.frequency = row_dict['duplication_frequency']
                interval.standard_error = row_dict['duplication_standard_error']
                interval.type = row_dict['type']
                interval.study = row_dict['study']
                interval.distance = interval.get_distance(
                    interval, Interval_base(chrom, start, end))
                result.append(interval)
        return result
