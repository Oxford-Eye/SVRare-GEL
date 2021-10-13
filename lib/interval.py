# given intervals, produce groups of intervals, where each member of a group overlaps with each other.
import numpy as np
import copy
from scipy.spatial import distance
from sklearn.cluster import DBSCAN


class Interval_base(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = end - start
        self.groups = []
        self._id = f'{chrom}-{start}-{end}'

    @staticmethod
    def overlap(this, other):
        if this.chrom == other.chrom:
            overlap_size = max(this.end, other.end) - \
                min(this.start, other.start)
            total_size = this.size + other.size
            if total_size > overlap_size:
                return True
        return False

    @staticmethod
    def get_distance(this, other):
        # 1 - (shared size / merged size)
        if this.chrom == other.chrom:
            distance = 1 - (min(this.end, other.end) - max(this.start, other.start)) / \
                (max(this.end, other.end) - min(this.start, other.start))

            return distance
        return None

    def __eq__(self, other):
        return ((self.chrom, self.start, self.end) == (other.chrom, other.start, other.end))

    def __lt__(self, other):
        if self.chrom < other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start < other.start:
                return True
            if self.start == other.start:
                if self.end < other.end:
                    return True
        return False

    def __le__(self, other):
        if self.chrom < other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start < other.start:
                return True
            if self.start == other.start:
                if self.end <= other.end:
                    return True
        return False

    def __gt__(self, other):
        if self.chrom > other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start > other.start:
                return True
            if self.start == other.start:
                if self.end > other.end:
                    return True
        return False

    def __ge__(self, other):
        if self.chrom > other.chrom:
            return True
        if self.chrom == other.chrom:
            if self.start > other.start:
                return True
            if self.start == other.start:
                if self.end >= other.end:
                    return True
        return False

    def __repr__(self):
        return self._id

    def __str__(self):
        return self._id

    @staticmethod
    def group(intervals):
        groups = []
        sorted_intervals = sorted(intervals)
        from datetime import datetime
        now = datetime.now().strftime("%H:%M:%S")
        print(f"{now} sorted")
        new_group = Group()
        for interval in sorted_intervals:
            check_interval = new_group.does_interval_belong_to_group(interval)
            if check_interval != False:
                new_group.add_interval(interval)
            else:
                groups.append(new_group)
                new_group = Group()
                new_group.add_interval(interval)
        groups.append(new_group)
        group_n = 0
        for interval in sorted_intervals:
            while group_n < len(groups) and interval in groups[group_n].intervals:
                group_n += 1
            while group_n < len(groups) and groups[group_n].does_interval_belong_to_group(interval):
                groups[group_n].add_interval(interval)
                group_n += 1
        return groups


class Group(object):
    def __init__(self, intervals=None):
        if intervals is None:
            self.intervals = []
        else:
            self.intervals = intervals
        self.chrom = None
        self.outer_start = None
        self.outer_end = None
        self.outer_id = None
        self.inner_start = None
        self.inner_end = None
        self.inner_id = None
        self.inner_interval = None
        self.outer_interval = None
        self.mean_start = None
        self.mean_end = None
        self.mean_interval = None
        self.mean_id = None
        self.stdev_start = None
        self.stdev_end = None

    def does_interval_belong_to_group(self, interval):
        # None: intervals empty, True: yes, False: no
        if not self.intervals:
            return None

        return Interval_base.overlap(interval, Interval_base(self.chrom, self.inner_start, self.inner_end))

    def add_interval(self, interval):
        #if self.intervals and not self.does_interval_belong_to_group(interval):
        #    print(
        #        f"interval: {interval} do not belong to group:{self} ")
        if self.intervals:
            self.outer_end = max(self.outer_end, interval.end)
            self.outer_start = min(self.outer_start, interval.start)
            self.inner_start = max(
                self.inner_start, interval.start)
            self.inner_end = min(self.inner_end, interval.end)
            self.outer_interval = Interval_base(
                self.chrom, self.outer_start, self.outer_end)
            self.inner_interval = Interval_base(
                self.chrom, self.inner_start, self.inner_end)
            self.mean_start = int(round(np.mean(
                [i.start for i in self.intervals] + [interval.start])))
            self.stdev_start = np.std(
                [i.start for i in self.intervals] + [interval.start])
            self.mean_end = int(round(np.mean(
                [i.end for i in self.intervals] + [interval.end])))
            self.stdev_end = np.std(
                [i.end for i in self.intervals] + [interval.end])
            self.mean_interval = Interval_base(
                self.chrom, self.mean_start, self.mean_end)
        else:
            self.chrom = interval.chrom
            self.outer_end = self.inner_end = self.mean_end = interval.end
            self.outer_start = self.inner_start = self.mean_start = interval.start
            self.outer_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.inner_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.mean_interval = Interval_base(
                interval.chrom, interval.start, interval.end)
            self.stdev_start = self.stdev_end = 0
        self.inner_id = f'{self.chrom}-{self.inner_start}-{self.inner_end}'
        self.mean_id = f'{self.chrom}-{self.mean_start}-{self.mean_end}'
        self.intervals.append(interval)
        self.intervals.sort()

    def produce_distance_matrix(self, method=None):
        # use start and size to produce distance matrix
        if method is None:
            method = interval_distance_for_cdist
        A = [(i.start, i.size) for i in self.intervals]

        return distance.cdist(A, A, method)

    def produce_distance_matrix_notyetimplemented(self, method=None):
        # use start and size to produce distance matrix
        from sklearn.neighbors import radius_neighbors_graph
        if method is None:
            method = interval_distance_for_cdist
        A = [(i.start, i.size) for i in self.intervals]

        return distance.cdist(A, A, method)

    def mean_distance_to_inner_interval(self):
        distances = [Interval_base.get_distance(self.inner_interval, i)
                     for i in self.intervals]
        return {'mean': np.mean(distances), 'stdev': np.std(distances)}

    def mean_distance_to_mean_interval(self):
        distances = [Interval_base.get_distance(self.mean_interval, i)
                     for i in self.intervals]
        return {'mean': np.mean(distances), 'stdev': np.std(distances)}

    def cluster(self, eps=0.1, min_samples=2):
        # making sub-groups with DBSCAN
        groups = {}
        if len(self.intervals) < 2:
            # no need to cluster
            group = Group()
            group.add_interval(self.intervals[0])
            groups[group.inner_id] = group
            return groups
        #print('compute D')
        #D = self.produce_distance_matrix(
        #    method=interval_distance_for_cdist)
        #print('D computed')
        #import sys
        #print(sys.getsizeof(D))
        #C = DBSCAN(eps=0.5, min_samples=2, metric='precomputed').fit(D)
        A = list(set([(i.start, i.size) for i in self.intervals]))
        '''
        with open('list', 'wt') as outf:
            for i in A:
                outf.write('\t'.join([str(j) for j in i]) + '\n')
        print(len(A))
        '''
        C = DBSCAN(eps=eps, metric=interval_distance_for_cdist, min_samples=min_samples).fit(A)

        for A_ind,A_val in enumerate(A):
            if C.labels_[A_ind] == -1:
                group = Group()
                for interval in (i for i in self.intervals if (i.start, i.size) == A_val):
                    group.add_interval(interval)
                groups[group.inner_id] = group
            else:
                if C.labels_[A_ind] not in groups:
                    groups[C.labels_[A_ind]] = Group()
                for interval in (i for i in self.intervals if (i.start, i.size) == A_val):
                    groups[C.labels_[A_ind]].add_interval(interval)


        return groups

    def normalise(self):
        pass

    def __eq__(self, other):
        return (self.inner_interval == other.inner_interval)

    def __lt__(self, other):
        return (self.inner_interval < other.inner_interval)

    def __le__(self, other):
        return (self.inner_interval <= other.inner_interval)

    def __gt__(self, other):
        return (self.inner_interval > other.inner_interval)

    def __ge__(self, other):
        return (self.inner_interval >= other.inner_interval)

    def __repr__(self):
        result = f'loose interval: {self.outer_interval}\n'
        result += f'core interval: {self.inner_interval}\n'
        result += f'mean interval: {self.mean_interval}\n'
        result += f'stdev start: {self.stdev_start}\n'
        result += f'stdev end: {self.stdev_end}\n'
        result += f'intervals: {self.intervals}\n'
        result += f''
        return result


def interval_distance_for_cdist_1(u, v):
    # customised distance measure
    # 1 - (shared size / merged size)
    distance = 1 - (min(u.start+u.size, v.start+v.size) - max(u.start, v.start)) / \
        (max(u.start+u.size, v.start+v.size) - min(u.start, v.start))

    return distance

def interval_distance_for_cdist(u, v):
    # customised distance measure
    # 1 - (shared size / merged size)
    distance = 1 - (min(u[0]+u[1], v[0]+v[1]) - max(u[0], v[0])) / \
        (max(u[0]+u[1], v[0]+v[1]) - min(u[0], v[0]))

    return distance
