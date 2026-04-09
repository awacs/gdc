# Convert a VCF/BCF file to eigenstrat format.
# Removes multiallelic and indel sites.
# usage: python vcf2eigenstrat.py -v vcf_file.vcf(.gz/.bcf) -o out_root
# will generate out_root.[snp,ind,geno].
# -i option is a .ind file to get population names and sex.
# --phased outputs phased 00/01/10/11 genotypes instead of standard 0/1/2.
# --aa polarizes by ancestral allele from the INFO AA field.

import argparse
import pysam

################################################################################

def parse_options():
    parser = argparse.ArgumentParser(
        description="Convert VCF/BCF to eigenstrat format (ind, snp, geno files)."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF/BCF file (.gz and .bcf supported)")
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

def get_aa(info):
    """
    Extract the ancestral allele from a pysam info dict.
    Handles both bare 'A' and 1000 Genomes 'A|||' pipe format.
    Returns a single uppercase nucleotide character, or None if absent/ambiguous.
    """
    try:
        aa = info["AA"]
        if isinstance(aa, (list, tuple)):
            aa = aa[0]
        aa = str(aa).split("|")[0].strip().upper()
        if aa in ("A", "T", "C", "G"):
            return aa
    except KeyError:
        pass
    return None

################################################################################

def decode_gt(gt_tuple, phased=False, flip=False):
    """
    Decode a pysam GT tuple into eigenstrat or phased output.
    gt_tuple elements are int allele indices (0=REF, 1=ALT) or None for missing.

    012 mode (default):
      Diploid: (0,0)->2, (0,1) or (1,0)->1, (1,1)->0, missing->9
      Haploid: (0,)->2, (1,)->0, missing->9
      flip inverts: 2<->0, 1 stays 1

    phased mode (--phased):
      Diploid: each allele written separately -> "00","01","10","11","99"
      Haploid: single character "0","1","9"
      flip swaps 0<->1 in each position
    """
    if gt_tuple is None:
        return "99" if phased else "9"

    def flip_allele(a):
        if a == 0: return 1
        if a == 1: return 0
        return a

    # Normalise to 0, 1, or None (treat any allele index >1 as alt)
    def norm(a):
        if a is None: return None
        return 0 if a == 0 else 1

    gt = tuple(norm(a) for a in gt_tuple)
    if flip:
        gt = tuple(flip_allele(a) if a is not None else None for a in gt)

    if len(gt) == 1:
        a = gt[0]
        if phased:
            return "9" if a is None else str(a)
        else:
            if a is None: return "9"
            return "2" if a == 0 else "0"

    if len(gt) == 2:
        a0, a1 = gt
        if phased:
            c0 = "9" if a0 is None else str(a0)
            c1 = "9" if a1 is None else str(a1)
            return c0 + c1
        else:
            if a0 is None or a1 is None:
                return "9"
            count = (1 if a0 == 0 else 0) + (1 if a1 == 0 else 0)
            return str(count)

    # Unexpected ploidy — treat as missing
    return "99" if phased else "9"

################################################################################

def main(options):
    """
    Convert VCF/BCF to eigenstrat format (ind, snp and geno files).
    """
    vcf = pysam.VariantFile(options.vcf)
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

    inds = list(vcf.header.samples)
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

    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos + 1  # pysam is 0-based; VCF POS is 1-based
        name = rec.id if rec.id else chrom + ":" + str(pos)
        ref = rec.ref
        alts = rec.alts or ()

        if len(alts) > 1:
            removed["multiallelic"] += 1
            continue

        alt = alts[0] if alts else "."

        if len(ref) != 1 or len(alt) != 1:
            removed["indel"] += 1
            continue

        # Ancestral allele polarization
        flip = False
        if options.aa:
            aa = get_aa(rec.info)
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
        snp.write("\t".join([name, chrom, "0.0", str(pos), snp_ref, snp_alt]) + "\n")

        # Write .geno
        geno_string = ""
        if options.ref:
            # Reference individual is hom-ref (or hom-ancestral after flip)
            geno_string = "00" if options.phased else "2"
        for indi in inds:
            gt_tuple = rec.samples[indi]["GT"]
            geno_string += decode_gt(gt_tuple, phased=options.phased, flip=flip)
        geno.write(geno_string + "\n")
        count += 1

    for f in [ind, snp, geno]:
        f.close()

    print(f"Done. Wrote {count} sites.")
    print(f"Excluded {sum(removed.values())} sites total.")
    for key, val in removed.items():
        if val:
            print(f"  Excluded {val} {key}")

    if options.aa and options.aa_unknown == "skip" and removed.get("unknownAA", 0) == count + sum(removed.values()):
        raise RuntimeError("--aa specified but no sites had an AA INFO field. All sites dropped. Check your VCF or use --aa-unknown ref.")

################################################################################

if __name__ == "__main__":
    options = parse_options()
    main(options)
