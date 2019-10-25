import numpy as np
from intervaltree import Interval

__version__ = "0.0.1"


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
            return Interval(int(start), int(end), chrom)
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
    return (a + a.T) / 2


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
    sim_mat = symmetrize(np.zeros((n, n)))
    print(sim_mat)
    # initialize `bed` objects
    beds = [Bed(b) for b in bedfiles]
    return beds
