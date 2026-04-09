# Converts a 4 haplotype macs output file (from msformatter)
# to msmc input format. Optionally add exponentially distributed phasing
# switch and flip errors between samples (0,1) and (2,3).

import sys
import argparse
import gdc
import numpy as np

################################################################################

def parse_options():
    parser = argparse.ArgumentParser(
        description="Convert ms/macs output to msmc/psmc input format."
    )
    parser.add_argument("-i", "--ms", required=True, help="Input ms/macs file (.gz supported)")
    parser.add_argument("-o", "--out", required=True, help="Output root")
    parser.add_argument("-c", "--chr", default="0", help="Chromosome name/number (default: 0)")
    parser.add_argument("-l", "--length", type=int, default=None, help="Sequence length (required unless --macs)")
    parser.add_argument("-s", "--switch-rate", type=float, default=0, dest="switch_rate",
                        help="Per-base phasing switch error rate")
    parser.add_argument("-f", "--flip-rate", type=float, default=0, dest="flip_rate",
                        help="Per-base phasing flip error rate")
    parser.add_argument("-m", "--msmc", action="store_true", help="Output msmc format")
    parser.add_argument("-p", "--psmc", action="store_true", help="Output psmc format")
    parser.add_argument("-a", "--macs", action="store_true", help="Input is macs format (infers length from header)")
    options = parser.parse_args()
    if not options.length and not options.macs:
        parser.error("Must specify --length unless using --macs input")
    return options

################################################################################

def add_phasing_errors(haps, switch_rate, flip_rate):
    """
    Add phasing errors to 2 haplotypes with switch and flip rates per site.
    """
    npos, nhap = haps.shape
    if nhap != 2:
        raise Exception("Can only add switch errors to 2 haplotypes")

    switches = np.random.uniform(size=npos) < switch_rate
    flips = np.random.uniform(size=npos) < flip_rate
    for i in range(npos):
        if switches[i]:
            tmp = haps[i:, 0].copy()
            haps[i:, 0] = haps[i:, 1]
            haps[i:, 1] = tmp
        if flips[i]:
            tmp = haps[i, 0]
            haps[i, 0] = haps[i, 1]
            haps[i, 1] = tmp

    return haps

################################################################################

def read_ms(ms_file, options):
    """
    Read a ms/macs file and return positions and haplotypes.
    """
    ms = gdc.open2(ms_file)

    nhap = length = None
    line = next(ms)
    if options.macs:
        nhap, length = [int(x) for x in line.split()[1:3]]
    else:
        length = options.length

    for line in ms:
        if line.startswith("segsites:"):
            npos = int(line.split()[1])
        elif line.startswith("positions:"):
            pos = np.array([int(length * float(p)) for p in line.split()[1:]])
            if len(pos) != npos:
                raise Exception("Number of positions does not match segsites")
            break

    haps = np.genfromtxt(ms, dtype=int, delimiter=1)   # rest of file is haplotypes
    haps = np.transpose(haps)

    if haps.shape[0] != npos:
        raise Exception("Number of positions doesn't match")
    if nhap and haps.shape[1] != nhap:
        raise Exception("Number of haplotypes doesn't match")

    return length, pos, haps

################################################################################

def main(options):
    length, pos, haps = read_ms(options.ms, options)

    npos, nhaps = haps.shape
    if nhaps not in [2, 4, 8]:
        raise Exception("Must have 2, 4 or 8 haplotypes")

    if options.switch_rate > 0 or options.flip_rate > 0:
        site_switch_rate = options.switch_rate * length / len(pos)
        site_flip_rate = options.flip_rate * length / len(pos)
        for i in range(nhaps // 2):
            haps[:, [2 * i, 2 * i + 1]] = add_phasing_errors(
                haps[:, [2 * i, 2 * i + 1]], site_switch_rate, site_flip_rate
            )

    alleles = np.zeros((npos, 2), dtype=str)
    alleles[:, 0] = "A"
    alleles[:, 1] = "T"

    if options.msmc:
        chrs = np.zeros(npos, dtype=str)
        chrs[:] = options.chr
        gdc.output_msmc(haps, chrs, pos, alleles, {"out": options.out})

    if options.psmc:
        gdc.output_psmc(haps, options.chr, pos, {"out": options.out})

################################################################################

if __name__ == "__main__":
    options = parse_options()
    main(options)
