import numpy as np
from interval import interval, inf, imath
from itertools import combinations

__version__ = "0.0.1"


class GenomicInterval(object):
    """
    Extend `interval` class to include contig/chromosome
    """

    def __init__(self, chrom, start, end):
        self.chr = chrom
        self.inf = int(start)
        self.sup = int(end)
        self.interval = interval([self.inf, self.sup])

    def __str__(self):
        return self.chr + ":" + str(self.inf) + "-" + str(self.sup)


class Bed:
    def __init__(self, fn=None):
        """
        Initialize object by storing file
        """
        self.file = fn
        self.fh = open(self.file, "r")

    def next(self):
        """
        Pop next interval from the BED file
        """
        newline = self.fh.readline()
        if newline != "":
            chrom, start, end = newline.rstrip().split("\t")
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
    return a + a.T - np.diag(np.diag(a))


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
        return min(
            [
                (a_int_b.sup - a_int_b.inf) / (a.sup - a.inf),
                (a_int_b.sup - a_int_b.inf) / (b.sup - b.inf),
            ]
        )
    else:
        return 0


def gis(bedfiles, names=None, prefix="similarity"):
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
    """
    # file handles for each BED file
    n = len(bedfiles)
    # store similarity matrix
    sim_mat = np.identity(n)
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
    print(sim_mat)
    print(minsim)
    return

