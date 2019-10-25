import os.path as path
import argparse


def validate(ARGS):
    """
    Validate command line arguments

    Parameters
    ----------
    ARGS : Namespace
        Command line arguments
    """
    # check for files existing
    nonex_files = [b for b in ARGS.bed if not path.exists(b)]
    if len(nonex_files) != 0:
        raise IOError(" ".join([nonex_files[0], "not found."]))
    # check `--names` and beds are of the same length
    if ARGS.names is not None:
        names = ARGS.names.split(",")
        if len(names) != len(ARGS.bed):
            raise IOError(
                "`--names` and `bed` must have the same number of arguments if `--names` is specified."
            )
    return {"bedfiles": ARGS.bed, "names": names, "prefix": ARGS.prefix}


def main():
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument("bed", type=str, help="BED file(s)", nargs="+")
    PARSER.add_argument(
        "-n", "--names", type=str, help="Comma-separated labels for BED files"
    )
    PARSER.add_argument(
        "-o",
        "--prefix",
        type=str,
        help="Prefix for output files.",
        default="similarity",
    )

    # parse arguments from command line
    ARGS = PARSER.parse_args()

    # validate command line arguments
    # import packages after parsing to speed up command line responsiveness
    validated_args = validate(ARGS)
    from .genomis import gis

    val = gis(**validated_args)
    if val is not None:
        print(val)
