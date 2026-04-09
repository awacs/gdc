# Convert a VCF file to eigenstrat format.
# Removes multiallelic and indel sites.
# usage: python vcf2eigenstrat.py -v vcf_file.vcf(.gz) -o out_root
# will generate out_root.[snp,ind,geno].
# -i option is a .ind file to get population names and sex.
# --phased outputs phased 00/01/10/11 genotypes instead of standard 0/1/2.
# --aa polarizes by ancestral allele from the INFO AA field.

import sys
import re
import argparse
import gdc

################################################################################

def parse_options():
    parser = argparse.ArgumentParser(
        description="Convert VCF to eigenstrat format (ind, snp, geno files)."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file (.gz supported)")
    parser.add_argument("-o", "--out", default="out", help="Output root (default: out)")
    parser.add_argument("-r", "--ref", default=None, help="Reference sample name (prepends a hom-ref row)")
    parser.add_argument("-i", "--ind", dest="indmap", default=None, help=".ind file for population/sex mapping")
    parser.add_argument("--indAsPop", action="store_true", help="Use individual name as population label")
    parser.add_argument("--phased", action="store_true", help="Output phased genotypes (00/01/10/11) instead of eigenstrat counts (0/1/2)")
    parser.add_argument("--aa", action="store_true", help="Polarize by ancestral allele from INFO AA field")
    parser.add_argument("--aa-unknown", choices=["skip", "ref"], default="skip",
                        help="Behaviour when AA is missing/ambiguous: 'skip' drops the site (default), 'ref' treats REF as ancestral")
    return parser.parse_args()

################################################################################

def get_aa(info_str):
    """
    Extract the ancestral allele from a VCF INFO string.
    Handles both bare 'AA=A' and 1000 Genomes 'AA=A|||' pipe format.
    Returns a single uppercase nucleotide character, or None if absent/ambiguous.
    """
    for field in info_str.split(";"):
        if field.startswith("AA="):
            aa = field[3:].split("|")[0].strip().upper()
            if aa in ("A", "T", "C", "G"):
                return aa
            return None
    return None

################################################################################

def decode_gt(gt_string, phased=False, flip=False):
    """
    Decode a VCF GT field into eigenstrat or phased output.

    012 mode (default):
      Diploid: 0/0->2, 0/1 or 1/0->1, 1/1->0, missing->9
      Haploid: 0->2, 1->0, missing->9
      flip inverts: 2<->0, 1 stays 1

    phased mode (--phased):
      Diploid: each allele written separately -> "00","01","10","11","99"
      Haploid: single character "0","1","9"
      flip swaps 0<->1 in each position
    """
    gt = gt_string.split(":")[0]
    alleles = re.split(r"[/|]", gt)

    def flip_allele(a):
        if a == "0":
            return "1"
        if a == "1":
            return "0"
        return a  # missing stays missing

    def to_char(a):
        """Normalise a single allele to '0', '1', or '.' for missing."""
        if a == "0":
            return "0"
        if a == "1":
            return "1"
        return "."

    alleles = [to_char(a) for a in alleles]

    if flip:
        alleles = [flip_allele(a) for a in alleles]

    if len(alleles) == 1:
        # Haploid
        a = alleles[0]
        if phased:
            return "9" if a == "." else a
        else:
            if a == ".":
                return "9"
            return "2" if a == "0" else "0"

    if len(alleles) == 2:
        a0, a1 = alleles
        if phased:
            c0 = "9" if a0 == "." else a0
            c1 = "9" if a1 == "." else a1
            return c0 + c1
        else:
            if a0 == "." or a1 == ".":
                return "9"
            count = (1 if a0 == "0" else 0) + (1 if a1 == "0" else 0)
            return str(count)

    # Unexpected ploidy — treat as missing
    return "99" if phased else "9"

################################################################################

def main(options):
    """
    Convert VCF to eigenstrat format (ind, snp and geno files).
    """
    vcf = gdc.open2(options.vcf)
    snp = open(options.out + ".snp", "w")
    ind = open(options.out + ".ind", "w")
    geno = open(options.out + ".geno", "w")

    removed = {"multiallelic": 0, "indel": 0}
    if options.aa:
        removed["unknownAA"] = 0
        removed["mismatchAA"] = 0
    count = 0

    pop_map = {}
    sex_map = {}
    if options.indmap:
        with open(options.indmap, "r") as ind_map_file:
            for line in ind_map_file:
                bits = line.rstrip("\n").split()
                pop_map[bits[0]] = bits[2]
                sex_map[bits[0]] = bits[1]

    for raw_line in vcf:
        # Handle bytes (gzipped input)
        line = raw_line.decode() if isinstance(raw_line, bytes) else raw_line

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            inds = line.split()[9:]
            if options.ref:
                ind.write(options.ref + "\tU\tREF\n")
            if options.indmap:
                for indi in inds:
                    ind.write(indi + "\t" + sex_map.get(indi, "U") + "\t" + pop_map.get(indi, "POP") + "\n")
            elif options.indAsPop:
                for indi in inds:
                    ind.write(indi + "\tU\t" + indi + "\n")
            else:
                for indi in inds:
                    ind.write(indi + "\tU\tPOP\n")
            continue

        bits = line.split()
        chrom, pos, name, ref, alt, info = bits[0], bits[1], bits[2], bits[3], bits[4], bits[7]

        if "," in alt:
            removed["multiallelic"] += 1
            continue
        if len(ref) != 1 or len(alt) != 1:
            removed["indel"] += 1
            continue

        if name == ".":
            name = chrom + ":" + pos

        # Ancestral allele polarization
        flip = False
        if options.aa:
            aa = get_aa(info)
            if aa is None:
                if options.aa_unknown == "skip":
                    removed["unknownAA"] += 1
                    continue
                # else: treat REF as ancestral, no flip
            elif aa == ref.upper():
                flip = False
            elif aa == alt.upper():
                flip = True
            else:
                removed["mismatchAA"] += 1
                continue

        # Write .snp (tab-separated): name, chrom, gpos, pos, allele1, allele2
        # allele1 is ancestral (or REF if not polarized); allele2 is derived (or ALT)
        snp_ref, snp_alt = (alt, ref) if flip else (ref, alt)
        snp.write("\t".join([name, chrom, "0.0", pos, snp_ref, snp_alt]) + "\n")

        # Write .geno
        geno_string = ""
        if options.ref:
            # Reference individual is hom-ref (or hom-ancestral after flip)
            geno_string = "00" if options.phased else "2"
        for gt in bits[9:]:
            geno_string += decode_gt(gt, phased=options.phased, flip=flip)
        geno.write(geno_string + "\n")
        count += 1

    for f in [ind, snp, geno]:
        f.close()

    print(f"Done. Wrote {count} sites.")
    print(f"Excluded {sum(removed.values())} sites total.")
    for key, val in removed.items():
        if val:
            print(f"  Excluded {val} {key}")

################################################################################

if __name__ == "__main__":
    options = parse_options()
    main(options)
