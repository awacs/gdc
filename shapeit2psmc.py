# Converts shapeit haps/sample files to psmc or msmc input files.
# For msmc, if you supply 2 individuals it will use both their haplotypes.
# If you supply 4, it will use one from each.

import sys
import argparse
import gzip
import gdc
import numpy as np

################################################################################

def parse_options():
    parser = argparse.ArgumentParser(
        description="Convert shapeit haps/sample files to psmc or msmc input."
    )
    parser.add_argument("-s", "--shapeit", required=True,
                        help="Root for .haps(.gz) and .sample files")
    parser.add_argument("-o", "--out", default="out", help="Output root (default: out)")
    parser.add_argument("-i", "--individuals", type=lambda s: s.split(","), default=[],
                        help="Comma-separated list of individuals from the sample file")
    parser.add_argument("-p", "--psmc", action="store_true", help="Output psmc format")
    parser.add_argument("-m", "--msmc", action="store_true", help="Output msmc format")
    return parser.parse_args()

################################################################################

def read_sample(sample_file, individuals):
    """
    Read the sample file and return the indices of the requested individuals.
    """
    names = []
    with open(sample_file) as f:
        next(f)
        next(f)                        # first two lines are header
        for line in f:
            names.append(line.split()[1])

    index = [names.index(i) for i in individuals]
    return index

################################################################################

def index_to_index(index):
    """
    Map individual indices to haplotype column indices.
    """
    if len(index) == 1:
        return np.array([2 * index[0], 2 * index[0] + 1])
    elif len(index) == 2:
        return np.array([2 * index[0], 2 * index[0] + 1, 2 * index[1], 2 * index[1] + 1])
    elif len(index) == 4:
        return np.array([2 * index[0], 2 * index[1], 2 * index[2], 2 * index[3]])
    else:
        raise Exception("Must be either 1, 2 or 4 individuals")

################################################################################

def main(options):
    if options.msmc and len(options.individuals) not in [2, 4]:
        raise Exception("msmc output needs either 2 or 4 individuals")
    if options.psmc and len(options.individuals) not in [1, 2]:
        raise Exception("psmc output needs either 1 or 2 individuals")

    index = read_sample(options.shapeit + ".sample", options.individuals)
    index = index_to_index(index)

    try:
        hap_file = gzip.open(options.shapeit + ".haps.gz", "rt")
    except IOError:
        hap_file = open(options.shapeit + ".haps", "r")

    chr = np.genfromtxt(hap_file, dtype=str, usecols=0)
    hap_file.seek(0)
    pos = np.genfromtxt(hap_file, dtype=int, usecols=2)
    hap_file.seek(0)
    alleles = np.genfromtxt(hap_file, dtype=str, usecols=(3, 4))
    hap_file.seek(0)
    haps = np.genfromtxt(hap_file, dtype=int, usecols=5 + index)

    if options.msmc:
        gdc.output_msmc(haps, chr, pos, alleles, {"out": options.out})

    if options.psmc:
        gdc.output_psmc(haps, chr[0], pos, {"out": options.out})

################################################################################

if __name__ == "__main__":
    options = parse_options()
    main(options)
