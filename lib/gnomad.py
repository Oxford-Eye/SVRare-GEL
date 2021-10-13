'''
parse gnomad
return an Interval instance, and a dictionary for later annotation
Note that the file has all sorts of SVs.
'''
import re
import gzip
import pysam
from lib import Interval_base


class Gnomad(object):
    tbx_header = [
        'chrom',
        'start',
        'end',
        'ID',
        'ref',
        'alt',
        'qual',
        'filter',
        'info',
    ]
    annotation_header = [
        'id',
        'cnv_type',
        'distance',
        'filter',
        'AN',
        'AC',
        'AF',
        'MALE_AN',
        'MALE_AC',
        'MALE_AF',
        'FEMALE_AN',
        'FEMALE_AC',
        'FEMALE_AF',
        'AFR_AN',
        'AFR_AC',
        'AFR_AF',
        'AFR_MALE_AN',
        'AFR_MALE_AC',
        'AFR_MALE_AF',
        'AFR_FEMALE_AN',
        'AFR_FEMALE_AC',
        'AFR_FEMALE_AF',
        'AMR_AN',
        'AMR_AC',
        'AMR_AF',
        'AMR_MALE_AN',
        'AMR_MALE_AC',
        'AMR_MALE_AF',
        'AMR_FEMALE_AN',
        'AMR_FEMALE_AC',
        'AMR_FEMALE_AF',
        'EAS_AN',
        'EAS_AC',
        'EAS_AF',
        'EAS_MALE_AN',
        'EAS_MALE_AC',
        'EAS_MALE_AF',
        'EAS_FEMALE_AN',
        'EAS_FEMALE_AC',
        'EAS_FEMALE_AF',
        'EUR_AN',
        'EUR_AC',
        'EUR_AF',
        'EUR_MALE_AN',
        'EUR_MALE_AC',
        'EUR_MALE_AF',
        'EUR_FEMALE_AN',
        'EUR_FEMALE_AC',
        'EUR_FEMALE_AF',
        'OTH_AN',
        'OTH_AC',
        'OTH_AF',
        'OTH_MALE_AN',
        'OTH_MALE_AC',
        'OTH_MALE_AF',
        'OTH_FEMALE_AN',
        'OTH_FEMALE_AC',
        'OTH_FEMALE_AF',
        'SVLEN',
        'ALGORITHMS',
        'EVIDENCE',
    ]

    pops = [
        'AFR',
        'AMR',
        'EUR',
        'OTH',
        'EAS',
    ]

    def __init__(self, fname, cnv_type):
        self.fname = fname
        self.tbx = pysam.TabixFile(fname)
        self.cnv_type = cnv_type
        # get info types
        # note that since the file doesn't contain genotype information
        # Number should only be one of ('.', '0', '1', 'A')
        allowed_info_types = {'.', '0', '1', 'A'}
        info_types = {}
        with gzip.open(fname, 'rt') as inf:
            for line in inf:
                if not line.startswith('##'):
                    break
                if line.startswith('##INFO'):
                    fields = line.split('<', 1)[1].split(',')
                    key = None
                    val = None
                    for field in fields:
                        if '=' in field:
                            field_key, field_val = field.split('=')
                            if field_key == 'ID':
                                key = field_val
                            elif field_key == 'Number':
                                if field_val not in allowed_info_types:
                                    raise ValueError(
                                        f"invalid Number:{field_val} for info field {field_key}")
                                val = field_val
                    if key is not None and val is not None:
                        info_types[key] = val
        self.info_types = info_types

    def get_gnomad(self, chrom, start, end):
        result = []
        try:
            it = self.tbx.fetch(chrom, max(0, start-1), end)
        except ValueError:
            return result
        for line in it:
            row_dict = dict(zip(self.tbx_header, line.split('\t')))
            info = self.parse_info(row_dict['info'])
            if info['SVTYPE'] in {'INS', 'BND', 'CPX'}:
                continue
            if self.cnv_type == 'LOSS' and info['SVTYPE'] == 'DUP':
                continue
            if self.cnv_type == 'GAIN' and info['SVTYPE'] == 'DEL':
                continue
            interval = Interval_base(
                row_dict['chrom'],
                int(row_dict['start']),
                int(info['END'])
            )
            interval.distance = interval.get_distance(
                interval, Interval_base(chrom, start, end))
            interval.cnv_type = self.cnv_type
            interval.filter = row_dict['filter']
            interval.annotation_header = self.annotation_header
            # set values for allelic unspecific fields, such as AN
            for key, val in info.items():
                if self.info_types[key] != 'A':
                    setattr(interval, key, val)
            if info['SVTYPE'] != 'MCNV':
                # not multiallelic
                for key, val in info.items():
                    if self.info_types[key] == 'A':
                        setattr(interval, key, val)
            else:
                # multiallelic
                for key in info:
                    if self.info_types[key] == 'A':
                        info[key] = info[key].split(',')
                AC = {
                    'AC': 0,
                    'FEMALE_AC': 0,
                    'MALE_AC': 0,
                }
                for pop in self.pops:
                    for suffix in ('AC', 'FEMALE_AC', 'MALE_AC'):
                        AC[f"{pop}_{suffix}"] = 0

                alts = row_dict['alt'].split(',')
                for ind in range(len(alts)):
                    cnv_number = int(
                        re.search(r'^<CN=(\d+)>$', alts[ind]).group(1))
                    if self.cnv_type == 'LOSS':
                        # LOSS
                        if (chrom == 'X' and cnv_number == 0):
                            AC['FEMALE_AC'] += int(info['FEMALE_AC'][ind])
                            AC['MALE_AC'] += int(info['MALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'FEMALE_AC', 'MALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                        elif chrom == 'X' and cnv_number == 1:
                            AC['FEMALE_AC'] += int(info['FEMALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'FEMALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                        elif chrom == 'Y' and cnv_number == 0:
                            AC['MALE_AC'] += int(info['MALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'MALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                        elif cnv_number < 2 and chrom not in ('X', 'Y'):
                            AC['FEMALE_AC'] += int(info['FEMALE_AC'][ind])
                            AC['MALE_AC'] += int(info['MALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'FEMALE_AC', 'MALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                    elif self.cnv_type == 'GAIN':
                        if chrom == 'X' and cnv_number == 2:
                            AC['MALE_AC'] += int(info['MALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'MALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                        elif chrom == 'Y' and cnv_number == 1:
                            AC['FEMALE_AC'] += int(info['FEMALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'FEMALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                        elif (chrom == 'Y' and cnv_number > 1) or cnv_number > 2:
                            AC['FEMALE_AC'] += int(info['FEMALE_AC'][ind])
                            AC['MALE_AC'] += int(info['MALE_AC'][ind])
                            for pop in self.pops:
                                for suffix in ('AC', 'FEMALE_AC', 'MALE_AC'):
                                    pop_key = f"{pop}_{suffix}"
                                    AC[pop_key] += int(info[pop_key][ind])
                    else:
                        raise ValueError(f"cnv type invalid: {self.cnv_type}")
                AC['AC'] = AC['FEMALE_AC'] + AC['MALE_AC']
                interval.AC = AC['AC']
                interval.FEMALE_AC = AC['FEMALE_AC']
                interval.MALE_AC = AC['MALE_AC']
                interval.AF = AC['AC'] / int(interval.AN)
                for gender in ('FEMALE', 'MALE'):
                    an = int(getattr(interval, f"{gender}_AN"))
                    if an == 0:
                        af = ''
                    else:
                        af = AC[f"{gender}_AC"] / an
                    setattr(interval, f"{gender}_AF", af)

                for pop in self.pops:
                    AC[f"{pop}_AC"] = AC[f"{pop}_FEMALE_AC"] + \
                        AC[f"{pop}_MALE_AC"]
                    for key in (f"{pop}_AC", f"{pop}_FEMALE_AC", f"{pop}_MALE_AC"):
                        setattr(interval, key, AC[key])
                    an = int(getattr(interval, f"{pop}_AN"))
                    if an == 0:
                        af = ''
                    else:
                        af = AC[f"{pop}_AC"] / an
                    setattr(interval, f"{pop}_AF", af)
                    for gender in ('MALE', 'FEMALE'):
                        an = int(getattr(interval, f"{pop}_{gender}_AN"))
                        if an == 0:
                            af = ''
                        else:
                            af = AC[f"{pop}_{gender}_AC"] / an
                        setattr(interval, f"{pop}_{gender}_AF", af)
            result.append(interval)
        return result

    def parse_info(self, info_raw):
        info = {}
        for info_field in info_raw.split(';'):
            if '=' in info_field:
                key, vals = info_field.split('=')
                info[key] = vals
            else:
                info[info_field] = True
        return info


def check_int(s):
    s = str(s)
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()
