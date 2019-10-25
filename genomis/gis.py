import numpy as np
from interval import interval
from itertools import combinations
import pandas as pd


class GenomicInterval(object):
    """
    Extend `interval` class to include contig/chromosome
    """

    def __init__(self, chrom, start, end):
        """
        Initialize object
        """
        self.chr = chrom
        self.inf = int(start)
        self.sup = int(end)
        self.size = self.sup - self.inf
        self.interval = interval([self.inf, self.sup])

    def __str__(self):
        """
        What to print when object is returned to interactive terminal
        """
        return self.chr + ":" + str(self.inf) + "-" + str(self.sup)

    def __repr__(self):
        """
        What to print when object is represented
        """
        return self.__str__()

    def __len__(self):
        """
        What to return when len() is called
        """
        return self.size


class Bed:
    def __init__(self, fn=None):
        """
        Initialize object by storing file
        """
        self.file = fn
        self.fh = open(self.file, "r")
        self.counter = -1

    def next(self):
        """
        Pop next interval from the BED file
        """
        newline = self.fh.readline()
        if newline != "":
            vals = newline.rstrip().split("\t")
            chrom = vals[0]
            start = vals[1]
            end = vals[2]
            self.counter += 1
            return GenomicInterval(chrom, start, end)
        return None

    def close(self):
        """
        Close object and release from memory
        """
        self.fh.close()


def symmetrize(a):
    """
    Symmetrize a matrix

    Parameters
    ----------
    a : numpy.array
        Matrix to be symmetrized
    """
    return a + a.T - np.diag(a.diagonal())


def similarity(a, b):
    """
    Calculate the symmetric overlap of two intervals

    Parameters
    ----------
    a : GenomicInterval
    b : GenomicInterval
    """
    if a.chr == b.chr:
        intersection = a.interval & b.interval
        a_int_b = GenomicInterval(a.chr, intersection[0].inf, intersection[0].sup)
        if len(a_int_b) == 0:
            return 0
        return min([len(a_int_b) / len(a), len(a_int_b) / len(b)])
    else:
        return 0


def gis(bedfiles, names=None, prefix="similarity", sim_thresh=0.5):
    """
    Calculate genomic similarity of BED files

    Parameters
    ----------
    beds : str
        Path to BED files to compare
    names : List<str> or NoneType
        Name for each of the input BED files
    prefix : str
        Output file prefix
    sim_thresh : float
        Minimum similarity threshold to consider recording a locus
    """
    # file handles for each BED file
    n = len(bedfiles)
    # store sample indexes as `names` is `names` is not defined
    if names is None:
        names = ["Sample_" + str(i) for i in range(n)]
    # store similarity matrix
    sim_mat = np.identity(n)
    # columns of information to store
    column_names = ["chr", "start", "end"] + names + ["similarity"]
    # initialize `bed` objects
    beds = [Bed(b) for b in bedfiles]
    # get first intervals from each file
    intvls = [b.next() for b in beds]
    # initialize the matrix
    for i, j in combinations(range(n), 2):
        sim_mat[i, j] = similarity(intvls[i], intvls[j])
    sim_mat = symmetrize(sim_mat)
    # store minimum similarity for each interval, and which sample it comes from
    minsim = [
        {"idx": j, "s": sim_mat[i, j]} for i, j in enumerate(np.argmin(sim_mat, axis=1))
    ]
    # records to keep for printing
    records = []
    # iterate over intervals from the sorted BED files
    while True:
        minsim_set = min([ms["s"] for ms in minsim])
        # record this set if the similarity of the set passes the threshold
        if minsim_set >= sim_thresh:
            chrom = intvls[0].chr
            # skip if not all intervals are on the same chromosome
            if np.all([intvl.chr == chrom for intvl in intvls]):
                # find interval that spans entire set of intervals
                hull = interval.hull([intvl.interval for intvl in intvls])
                set_locus = GenomicInterval(chrom, hull[0].inf, hull[0].sup)
                records.append(
                    dict(
                        (colname, v)
                        for colname, v in zip(
                            column_names,
                            [
                                set_locus.chr,
                                set_locus.inf,
                                set_locus.sup,
                                *[b.counter for b in beds],
                                minsim_set,
                            ],
                        )
                    )
                )
        # find sample with the smallest upper bound
        update_idx = np.argmin([intvl.sup for intvl in intvls])
        # pop this interval
        intvls[update_idx] = beds[update_idx].next()
        # check that we're not at the end of the file
        if intvls[update_idx] is None:
            break
        # recalculate column of sim_mat (calc once, ensure sim_mat is symmetric)
        for i in range(update_idx):
            sim_mat[i, update_idx] = similarity(intvls[i], intvls[update_idx])
            sim_mat[update_idx, i] = sim_mat[i, update_idx]
        for i in range(update_idx + 1, n):
            sim_mat[i, update_idx] = similarity(intvls[i], intvls[update_idx])
            sim_mat[update_idx, i] = sim_mat[i, update_idx]
        # update minsim for (update_idx)-th sample
        minsim[update_idx]["idx"] = np.argmin(sim_mat[:, update_idx])
        minsim[update_idx]["s"] = sim_mat[update_idx, minsim[update_idx]["idx"]]
        # update minsim for any sample where minsim[j]["idx"] == update_idx
        for i in [i for i, ms in enumerate(minsim) if ms["idx"] == update_idx]:
            minsim[i]["idx"] = np.argmin(sim_mat[:, i])
            minsim[i]["s"] = sim_mat[i, minsim[i]["idx"]]
    # save records as a DataFrame
    df = pd.DataFrame(records, columns=column_names)
    # save to output
    df.to_csv(prefix + ".tsv", index=False, sep="\t")
    return df
