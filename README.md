# gdc: Genetic Data Conversion

A collection of Python 3 scripts for converting genetic data between common population genetics file formats. All scripts require Python 3.6+.

## Shared Library

### `gdc.py`

Utility module imported by most scripts. Provides:

- **`open2(file, mode="r")`** — Opens a regular file or a gzipped file (detected by `.gz` extension). Gzipped files automatically open in text mode (`"rt"`) unless binary mode is explicitly requested. This is the standard way all scripts handle transparent `.gz` input.
- **`output_msmc(haps, chr, pos, alleles, options)`** — Writes an `.msmc` file from a numpy haplotype array (4 or 8 columns). When 8 haplotypes are present, uses columns `[0,2,4,6]` (one per individual). Output columns: chromosome, position, distance from last site, allele string.
- **`output_psmc(haps, chr, pos, options)`** — Writes a `.psmc` file from a numpy haplotype array (2 or 4 columns). Encodes heterozygosity in 100bp blocks as `T` (hom) or `K` (het), wrapped at 60 characters per line.

---

## Scripts

### `vcf2eigenstrat.py`

Converts a VCF file to eigenstrat format (`.snp`, `.ind`, `.geno` files) using standard allele-count encoding.

**Usage:**
```
python vcf2eigenstrat.py -v input.vcf.gz -o output_root [options]
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-v/--vcf` | Input VCF file (`.gz` supported). Required. |
| `-o/--out` | Output file root (default: `out`). Produces `{root}.snp`, `{root}.ind`, `{root}.geno`. |
| `-r/--ref` | Reference sample name. Prepends a homozygous-reference row (genotype `2`) to every line in `.geno`. |
| `-i/--ind` | `.ind` file mapping individuals to populations and sex. Three whitespace-delimited columns: name, sex, population. |
| `--indAsPop` | Use each individual's name as its population label. |

**Input:** Standard VCF (v4.x). Only biallelic SNPs are kept; multiallelic sites (comma in ALT) and indels (REF or ALT length > 1) are removed and counted.

**Output:**
- `.snp` — Tab-separated: SNP name, chromosome, genetic position (0.0), physical position, REF allele, ALT allele. If the VCF has no rsID (`.`), uses `chr:pos` as the name.
- `.ind` — Tab-separated: individual name, sex (`U` if unknown), population label.
- `.geno` — One row per SNP, one character per individual. Eigenstrat encoding: `2` = homozygous reference, `1` = heterozygous, `0` = homozygous alternate, `9` = missing. Haploid genotypes (e.g. chrX in males): `0` → `2`, `1` → `0`, missing → `9`.

**Logic:**
1. Parse the `#CHROM` header to extract individual names and write `.ind`.
2. For each data line, split on whitespace, check ALT for commas (multiallelic) and REF/ALT length (indels).
3. Extract the GT field (first colon-delimited element of each sample column), decode to allele count, and concatenate into the `.geno` line.

---

### `vcf2eigenstrat.py` — extended options

In addition to the base arguments above, the following flags are supported:
| Flag | Description |
|------|-------------|
| `--phased` | Output phased genotypes (`00`/`01`/`10`/`11`) instead of allele counts (`0`/`1`/`2`). Haploid sites output a single character. |
| `--aa` | Polarize by ancestral allele. Reads the `AA` tag from the VCF INFO field. When AA matches ALT, flips REF/ALT in the `.snp` file and inverts genotype encoding. |
| `--aa-unknown` | Behaviour when AA is missing or ambiguous: `skip` (default, drops the site) or `ref` (treats REF as ancestral, no flip). |

**Ancestral allele parsing (`get_aa`):**
Searches the INFO field (column 8) for a tag starting with `AA=`. Handles the 1000 Genomes pipe-delimited format (`AA=A|||`) by splitting on `|` and taking the first element. Returns `None` (triggering skip or ref-fallback) if the value is not a clear nucleotide (`A`, `T`, `C`, `G`).

**Flip logic:**
When `--aa` is active and AA matches the ALT allele:
- `.snp` columns swap: allele1 becomes ALT (now ancestral), allele2 becomes REF (now derived).
- In `--phased` mode: `0` and `1` swap in each allele position.
- In default (012) mode: `2` ↔ `0`, `1` stays `1`.
- Sites where AA matches neither REF nor ALT are counted as `mismatchAA` and skipped.

---

### `eigenstrat2vcf.py`

Converts eigenstrat format back to VCF. Writes to stdout.

**Usage:**
```
python eigenstrat2vcf.py -r data_root [options] > output.vcf
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-r/--root` | Root for eigenstrat files (`{root}.snp`, `{root}.geno`, `{root}.ind`). Required. |
| `-i/--inds` | File listing individual names to include (one per line). |
| `-p/--pops` | File listing populations to include (one per line). |
| `-s/--snps` | File listing SNP names to include (one per line). |

**Requires:** `pyEigenstrat` package.

**Input:** Standard eigenstrat triplet (`.snp`, `.ind`, `.geno`), packed or unpacked.

**Output:** VCF v4.0 to stdout. Genotypes are mapped: `2` → `0/0`, `1` → `0/1`, `0` → `1/1`, `9` → `./.`.

**Logic:**
1. Load data via `pyEigenstrat.load()`, optionally filtering by individuals, populations, or SNPs.
2. Write VCF header lines, then iterate over SNPs, converting each eigenstrat genotype row back to VCF GT fields.

---

### `arlequin2eigenstrat.py`

Converts an Arlequin `.arp` file (generated by fastsimcoal2 with `-s0`) to eigenstrat format.

**Usage:**
```
python arlequin2eigenstrat.py -a input.arp -o output_root [-p]
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-a/--arp` | Input `.arp` file (`.gz` supported). Required. |
| `-o/--out` | Output root (default: `out`). |
| `-p/--phased` | Keep output phased (one haplotype per column). Default behaviour unphases by summing pairs of haplotypes into diploid genotype counts. |

**Input:** Arlequin `.arp` format with DNA sequence data (`-s0` option in fastsimcoal2). Supports multiple populations and chromosomes. Automatically detects ascertained data sections.

**Output:** Eigenstrat triplet. SNP names are `{chrom}_{pos}`. Population labels come from the `SampleName` fields in the `.arp` file. In unphased mode (default), genotype values are `0`/`1`/`2`; in phased mode, `0`/`1`.

**Logic:**
1. **`load_from_arp`**: Parses the `.arp` file structure, reading chromosome site positions, sample names, sample sizes, and genotype strings. Genotype values `{1,2,3}` are collapsed to `1` (derived), `0` stays `0` (ancestral).
2. **`unphase`**: For each population, pairs consecutive haplotype columns and sums them to produce diploid counts (`0`/`1`/`2`). Halves the sample count.
3. Writes `.snp`, `.ind`, `.geno` from the parsed data using numpy's `savetxt`.

---

### `chromopainter2eigenstrat.py`

Converts a ChromoPainter-format haplotype file to eigenstrat format.

**Usage:**
```
python chromopainter2eigenstrat.py --hap hap_file -d donor_file -r recipient_file -o output_root
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `--hap` | Haplotype file (`.gz` supported). One site per row, one haplotype per column. Required. |
| `-d/--donor` | Donor population file (population name + haplotype count per line). Required. |
| `-r/--recipient` | Recipient population file (same format). Required. |
| `-o/--out` | Output root (default: `out`). |

**Input:**
- Haplotype file: columns are `chrom pos allele1 allele2 ...`. Columns are ordered donors first, then recipients.
- Donor/recipient files: each line is `population_name num_haplotypes`. Haplotype counts must be even (paired into diploid individuals).

**Output:** Eigenstrat triplet. Genotype is the count of the first observed allele per site. SNP names are `chr:pos`. Multiallelic sites (>2 alleles) are skipped.

**Logic:**
1. Read donor and recipient files to build the `.ind` file, pairing haplotypes into diploid individuals.
2. For each haplotype line, identify the two alleles present, sum matches to the first allele per individual pair, and write to `.geno`.

---

### `vcf2freq.py`

Extracts per-population allele frequencies from a VCF.

**Usage:**
```
python vcf2freq.py -i input.vcf -p panel_file > freqs.tsv
```
Reads from stdin if `-i` is omitted.

**Arguments:**
| Flag | Description |
|------|-------------|
| `-i/--input` | Input VCF (default: stdin). |
| `-p/--panel` | Two-column file mapping sample IDs to population labels. Required. |

**Input:** VCF with phased genotypes (expects `0|0`, `0|1`, etc. — reads character positions `[0]` and `[1]` of each GT field). Panel file: `sample_id population` per line.

**Output:** Tab-separated table to stdout. Header: `SNPID CHR POS REF ALT pop1 pop2 ...`. Each row is one SNP with ALT allele frequency per population (4 decimal places). Sites where any population has zero callable alleles are skipped.

**Logic:**
1. Read panel file into a sample→population map.
2. For each VCF data line, count ALT alleles and total callable alleles per population by reading character positions 0 and 1 of each sample's GT string.
3. Compute frequency = ALT count / total. Skip on `ZeroDivisionError`.

---

### `vcf2hetfa.py`

Converts a VCF to heterozygous FASTA (hetfa) format for PSMC input, or to two phased FASTA files.

**Usage:**
```
# Hetfa mode (default) — pipe to fq2psmcfa:
python vcf2hetfa.py -v input.vcf -r ref.fa -s sample_name -c chr1 | fold | fq2psmcfa - > sample.psmcfa

# Haplotype mode — two gzipped FASTA files:
python vcf2hetfa.py -v input.vcf -r ref.fa -s sample_name -c chr1 -h -o output_root
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-v/--vcf` | Input VCF. Required. |
| `-r/--ref` | Reference FASTA (indexed by pyfaidx). Required. |
| `-s/--sample` | Sample name (must match a VCF column). Required. |
| `-c/--chrom` | Chromosome to extract. Required. |
| `-o/--out` | Output root. Required for `-h` mode; if omitted in hetfa mode, writes to stdout. |
| `-h/--haplotypes` | Output two phased `.fa.gz` files instead of a single hetfa. |
| `-m/--mask` | Mask FASTA file (digits 0-9 per position). |
| `-a/--mask_value` | Minimum mask value to keep a site (must be specified with `-m`). |

**Requires:** `pyfaidx` package.

**Input:** VCF with phased genotypes (`0|0`, `0|1`, `1|0`, `1|1`) and an indexed reference FASTA.

**Output:**
- **Hetfa mode**: Single FASTA to stdout (or gzipped file). Homozygous ref → ref base, homozygous alt → alt base, heterozygous → IUPAC ambiguity code (`R`, `Y`, `S`, `W`, `K`, `M`). Non-biallelic/unphased/missing → `N`.
- **Haplotype mode**: Two gzipped FASTA files (`{root}.0.fa.gz`, `{root}.1.fa.gz`), one per phased chromosome. Each base is either REF or ALT at variant sites, with reference sequence filling in between.

**Masking:** When `-m` and `-a` are provided, any position where the mask value is below the threshold (or non-numeric) is replaced with `N` in the output.

**Logic:**
1. Open the reference FASTA and (optionally) the mask FASTA via pyfaidx.
2. Walk through VCF data lines in order. For each variant, emit the reference sequence from the last position up to this one, then emit the genotype-appropriate base(s).
3. After the last variant, emit the remaining reference tail.

---

### `ms2psmc.py`

Converts ms/macs coalescent simulation output to MSMC or PSMC input format.

**Usage:**
```
python ms2psmc.py -i sim_output.txt -o output_root -l 1000000 [-m] [-p] [options]
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-i/--ms` | Input ms/macs file (`.gz` supported). Required. |
| `-o/--out` | Output root. Required. |
| `-l/--length` | Sequence length in bp. Required unless `--macs`. |
| `-c/--chr` | Chromosome name (default: `0`). |
| `-m/--msmc` | Output MSMC format. |
| `-p/--psmc` | Output PSMC format. |
| `-a/--macs` | Input is macs format (reads length from header). |
| `-s/--switch-rate` | Per-base phasing switch error rate (default: 0). |
| `-f/--flip-rate` | Per-base phasing flip error rate (default: 0). |

**Input:** Standard ms/macs output: a header line, a `segsites:` line, a `positions:` line, then one row per haplotype (0/1 characters with no delimiter). Must contain 2, 4, or 8 haplotypes.

**Output:**
- **MSMC** (`.msmc`): Tab-separated: chromosome, position, distance from previous site, allele string (one character per haplotype). Uses columns `[0,2,4,6]` if 8 haplotypes are present.
- **PSMC** (`.psmc`): 100bp-block encoding where `T` = homozygous block, `K` = block containing a het site.

**Phasing errors:** When `--switch-rate` or `--flip-rate` > 0, errors are introduced between haplotype pairs (0,1), (2,3), etc. Switch errors swap all downstream alleles between the pair; flip errors swap a single site.

**Logic:**
1. **`read_ms`**: Parse the ms/macs file. Convert fractional positions to absolute bp using the sequence length. Read haplotype rows via `numpy.genfromtxt` with `delimiter=1`.
2. Optionally add phasing errors to each haplotype pair.
3. Delegate output to `gdc.output_msmc()` or `gdc.output_psmc()`.

---

### `shapeit2psmc.py`

Converts SHAPEIT2 phased output (`.haps`/`.sample` files) to MSMC or PSMC input format.

**Usage:**
```
python shapeit2psmc.py -s shapeit_root -i ind1,ind2 -o output_root [-m] [-p]
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-s/--shapeit` | Root for SHAPEIT files (`{root}.haps` or `{root}.haps.gz`, and `{root}.sample`). Required. |
| `-o/--out` | Output root (default: `out`). |
| `-i/--individuals` | Comma-separated list of individual IDs from the `.sample` file. |
| `-p/--psmc` | Output PSMC format (needs 1 or 2 individuals). |
| `-m/--msmc` | Output MSMC format (needs 2 or 4 individuals). |

**Input:**
- `.haps`: SHAPEIT haplotype format. Columns: chromosome, SNP ID, position, allele1, allele2, then one column per haplotype (0/1).
- `.sample`: SHAPEIT sample file. Two header lines, then one line per individual. The second column is the individual ID.

**Haplotype selection:**
- 1 individual → 2 haplotypes (both chromosomes).
- 2 individuals → 4 haplotypes (both chromosomes of each).
- 4 individuals → 4 haplotypes (one chromosome from each).

**Output:** MSMC (`.msmc`) or PSMC (`.psmc`) format via `gdc.output_msmc()` / `gdc.output_psmc()`.

**Logic:**
1. **`read_sample`**: Parse the `.sample` file to find column indices for the requested individuals.
2. **`index_to_index`**: Map individual indices to haplotype column indices in the `.haps` file (each individual has 2 consecutive columns starting at column 5).
3. Load the full `.haps` file with `numpy.genfromtxt` (4 passes: chr, pos, alleles, haplotypes), then delegate to gdc output functions.

---

### `maskfa.py`

Applies a quality mask to a FASTA file, replacing low-quality positions with `N`.

**Usage:**
```
python maskfa.py -f input.fa -m mask.fa -c 20 > masked.fa
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-f/--fasta` | Input FASTA file. Required. |
| `-m/--mask` | Mask FASTA file (digits 0-9 per position). If omitted, outputs a binary callability track. |
| `-c/--level` | Minimum quality level to keep (default: 1). |

**Requires:** `pyfaidx` package.

**Input:** A FASTA file and a corresponding mask FASTA where each position is a digit (0-9) representing quality/coverage.

**Output (to stdout):**
- **With mask**: The original FASTA with positions below the threshold replaced by `N`. Non-digit mask characters also trigger `N`.
- **Without mask**: A binary callability FASTA where each position is `1` if it's a called base (A/C/G/T or IUPAC het code) or `0` otherwise.

Output is wrapped at 50 characters per line.

---

### `polysites2vcf.py`

Converts Shop Mallick's polysite format (SGDP project) to VCF. Highly specific to that dataset's format.

**Usage:**
```
python polysites2vcf.py -i polysites_file [-c chr1] > output.vcf
```

**Arguments:**
| Flag | Description |
|------|-------------|
| `-i/--input` | Input polysite file (default: stdin). |
| `-c/--chrom` | Restrict output to a single chromosome. |

**Input:** SGDP polysite format. Header lines starting with `##` contain sample metadata. A `#` line marks the end of the header. Data lines contain chromosome, position, reference base, and IUPAC-coded genotypes for each population group.

**Output:** VCF v4.2 to stdout. IUPAC codes are decoded to diploid genotypes (e.g., `R` → `A/G` → `0/1`). Monomorphic sites are skipped.

---

## Dependencies

- **Python 3.6+**
- **NumPy** — used by `arlequin2eigenstrat.py`, `ms2psmc.py`, `shapeit2psmc.py`
- **pyfaidx** — used by `vcf2hetfa.py`, `maskfa.py`
- **pyEigenstrat** — used by `eigenstrat2vcf.py`

## File Format Quick Reference

| Format | Files | Description |
|--------|-------|-------------|
| Eigenstrat | `.snp`, `.ind`, `.geno` | SNP info, individual info, genotype matrix (0/1/2/9) |
| VCF | `.vcf` / `.vcf.gz` | Variant Call Format |
| MSMC | `.msmc` | Input for MSMC (multiple sequentially Markovian coalescent) |
| PSMC | `.psmc` | Input for PSMC (pairwise SMC) |
| Arlequin | `.arp` | Population genetics simulation output |
| SHAPEIT | `.haps`, `.sample` | Phased haplotype data |
| hetfa | `.hetfa.fa.gz` | Heterozygous FASTA for PSMC |
