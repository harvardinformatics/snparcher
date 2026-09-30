#!/usr/bin/env python3
"""Generate the simulated fixture for the RepAdapt equivalence tests.

Writes, under the output directory (default: this script's directory):

    reference.fasta              three contigs, about 24.5 kb in total
    genes.gff                    a few gene records (RepAdapt requires --gff_file)
    fastq/S{1..4}_{1,2}.fastq.gz paired 2x150 reads, at most 1,000 pairs per sample
    truth.tsv                    the simulated variants and genotypes

The generator is pure Python and seeded, and uses only ``Random.random()`` so
its output does not depend on the Python version. Rerunning it reproduces the
same sequences; gzip members are written with mtime 0.

What each feature is for:

- An exact 1 kb repeat in ctgA gives MAPQ-0 pairs, which RepAdapt's
  ``samtools view -q 10`` removes.
- A near-identical copy of a 2 kb ctgA segment in ctgB, with one difference
  every 500 bp, gives pairs that span one difference MAPQ ~27 in bwa. They
  pass ``-q 10`` but pull INFO/MQ below 30, so SNPs there fail ``MQ<30``.
- SNPs that are hom-alt in every sample fail ``AC=AN``.
- Indels of 1-30 bp, some inside short tandem repeats (where bwa places the
  gap inconsistently between reads), each with a SNP 8-20 bp away, give indel
  realignment work that changes the pileup at those SNPs.
- About 5% of pairs are PCR duplicates (a fragment sequenced again).
- About 3% of pairs have inserts shorter than the read, so the reads run into
  adapter sequence.
- About 2% of pairs have an unmappable mate. After the MAPQ filter, fixmate
  turns the mapped read into a single-end read.

Read names are SRA-like (``SIM.<n>``). Names starting with @NS, @NB or @A0
switch on fastp's polyG trimming, and Illumina tile/x/y names switch on
Picard's optical-duplicate model; both are avoided.

At most 1,000 pairs per sample keeps RepAdapt's multi-threaded fastp
deterministic: fastp 0.20.1 writes packs of 1,000 pairs in completion order.
"""

import argparse
import gzip
import math
from pathlib import Path
from random import Random

SEED = 20260930

CONTIG_LENGTHS = {"ctgA": 15000, "ctgB": 8000, "ctgC": 1500}

# Exact repeat: ctgA[EXACT_SRC] is copied over ctgA[EXACT_DEST].
EXACT_SRC = 2000
EXACT_DEST = 9000
EXACT_LEN = 1000

# Near-identical copy: ctgA[NEAR_SRC:NEAR_SRC + NEAR_LEN] is copied over
# ctgB[NEAR_DEST:...], then one base in each NEAR_DIFF_SPACING-bp block is
# changed so the two copies differ at NEAR_LEN / NEAR_DIFF_SPACING sites.
NEAR_SRC = 11500
NEAR_DEST = 3000
NEAR_LEN = 2000
NEAR_DIFF_SPACING = 500

SAMPLES = ["S1", "S2", "S3", "S4"]
PAIRS_PER_SAMPLE = 1000
READ_LEN = 150
INSERT_MEAN = 350
INSERT_SD = 50
INSERT_MIN = 200
INSERT_MAX = 500
SHORT_INSERT_MIN = 100
SHORT_INSERT_MAX = 140
DUPLICATE_FRACTION = 0.05
SHORT_INSERT_FRACTION = 0.03
JUNK_MATE_FRACTION = 0.02
ERROR_RATE = 0.005

# TruSeq adapters as they appear in R1 and R2 read-through.
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

GENES = [
    ("ctgA", 501, 1800, "+"),
    ("ctgA", 4001, 6500, "-"),
    ("ctgA", 12001, 13200, "+"),
    ("ctgB", 801, 2400, "+"),
    ("ctgB", 5501, 7200, "-"),
    ("ctgC", 201, 1200, "+"),
]


class Rng:
    """Seeded helpers built only on Random.random(), which is stable across
    Python versions (unlike randrange/shuffle/gauss/choices)."""

    def __init__(self, seed):
        self._rng = Random(seed)

    def random(self):
        return self._rng.random()

    def below(self, n):
        return min(int(self._rng.random() * n), n - 1)

    def between(self, lo, hi):
        """Integer in [lo, hi]."""
        return lo + self.below(hi - lo + 1)

    def base(self, exclude=None):
        choices = [b for b in BASES if b != exclude]
        return choices[self.below(len(choices))]

    def seq(self, n):
        return "".join(BASES[self.below(4)] for _ in range(n))

    def gauss(self, mu, sigma):
        u1 = 1.0 - self._rng.random()
        u2 = self._rng.random()
        return mu + sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

    def shuffle(self, items):
        for i in range(len(items) - 1, 0, -1):
            j = self.below(i + 1)
            items[i], items[j] = items[j], items[i]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def build_reference(rng):
    ref = {name: list(rng.seq(length)) for name, length in CONTIG_LENGTHS.items()}

    ctga = ref["ctgA"]
    ctga[EXACT_DEST : EXACT_DEST + EXACT_LEN] = ctga[EXACT_SRC : EXACT_SRC + EXACT_LEN]

    near = ctga[NEAR_SRC : NEAR_SRC + NEAR_LEN]
    near_diffs = []
    for block_start in range(0, NEAR_LEN, NEAR_DIFF_SPACING):
        offset = block_start + NEAR_DIFF_SPACING // 2 + rng.between(-50, 50)
        near[offset] = rng.base(exclude=near[offset])
        near_diffs.append(offset)
    ref["ctgB"][NEAR_DEST : NEAR_DEST + NEAR_LEN] = near

    return {name: "".join(bases) for name, bases in ref.items()}, near_diffs


def repeat_intervals():
    """0-based half-open intervals where only deliberate variants may go."""
    return [
        ("ctgA", EXACT_SRC, EXACT_SRC + EXACT_LEN),
        ("ctgA", EXACT_DEST, EXACT_DEST + EXACT_LEN),
        ("ctgA", NEAR_SRC, NEAR_SRC + NEAR_LEN),
        ("ctgB", NEAR_DEST, NEAR_DEST + NEAR_LEN),
    ]


def in_intervals(contig, pos, pad, intervals):
    return any(c == contig and s - pad <= pos < e + pad for c, s, e in intervals)


def random_genotypes(rng, allow_all_alt=False):
    """Diploid genotypes for every sample, with at least one alt allele and
    (unless allowed) at least one ref allele across the cohort."""
    while True:
        gts = []
        for _ in SAMPLES:
            r = rng.random()
            if r < 0.4:
                gts.append((0, 0))
            elif r < 0.8:
                gts.append((0, 1) if rng.random() < 0.5 else (1, 0))
            else:
                gts.append((1, 1))
        alt = sum(a + b for a, b in gts)
        if alt == 0:
            continue
        if not allow_all_alt and alt == 2 * len(SAMPLES):
            continue
        return gts


def het_genotypes(rng):
    """Every sample het or hom-ref, at least two het: gives a well-supported
    site whose only failing filter is MQ."""
    while True:
        gts = []
        for _ in SAMPLES:
            if rng.random() < 0.7:
                gts.append((0, 1) if rng.random() < 0.5 else (1, 0))
            else:
                gts.append((0, 0))
        if sum(a + b for a, b in gts) >= 2:
            return gts


def make_variant(rng, ref, contig, pos, kind, size=1, genotypes=None, tag="snp", inserted=None):
    ref_base = ref[contig][pos]
    if kind == "snp":
        ref_allele = ref_base
        alt_allele = rng.base(exclude=ref_base)
    elif kind == "del":
        ref_allele = ref[contig][pos : pos + size + 1]
        alt_allele = ref_base
    elif kind == "ins":
        ref_allele = ref_base
        alt_allele = ref_base + (inserted if inserted is not None else rng.seq(size))
    else:
        raise ValueError(kind)
    return {
        "contig": contig,
        "pos": pos,
        "ref": ref_allele,
        "alt": alt_allele,
        "kind": kind,
        "tag": tag,
        "gts": genotypes if genotypes is not None else random_genotypes(rng),
    }


def place_variants(rng, ref, near_diffs):
    variants = []
    occupied = []  # (contig, start, end) with padding already applied
    blocked = repeat_intervals()

    def free(contig, start, end, pad):
        if start < 200 or end > CONTIG_LENGTHS[contig] - 200:
            return False
        return not any(c == contig and s < end + pad and start < e + pad for c, s, e in occupied)

    def claim(contig, start, end):
        occupied.append((contig, start, end))

    # Indels, each with a SNP 8-20 bp downstream that shares its carriers, so
    # realignment changes the pileup at the SNP: short and long indels in
    # unique sequence, and indels of one repeat unit inside short tandem
    # repeats, where bwa places the gap inconsistently between reads. The
    # tandem repeats are planted into the reference here.
    indel_specs = [
        ("del", 1, None),
        ("ins", 2, None),
        ("del", 3, None),
        ("ins", 5, None),
        ("del", 8, None),
        ("ins", 10, None),
        ("del", 15, None),
        ("ins", 18, None),
        ("del", 25, None),
        ("ins", 30, None),
        ("del", 1, "A" * 10),
        ("ins", 1, "T" * 9),
        ("del", 2, "CA" * 6),
        ("ins", 2, "GT" * 7),
        ("del", 3, "TTG" * 5),
        ("ins", 4, "ACAG" * 4),
    ]
    for kind, size, repeat in indel_specs:
        span = len(repeat) if repeat else size
        while True:
            contig = ["ctgA", "ctgA", "ctgA", "ctgB", "ctgB", "ctgC"][rng.below(6)]
            pos = rng.between(200, CONTIG_LENGTHS[contig] - 280)
            snp_pos = pos + 1 + span + rng.between(8, 20)
            if in_intervals(contig, pos, 400, blocked):
                continue
            if not free(contig, pos, snp_pos + 1, 40):
                continue
            break
        inserted = None
        tag = f"indel_{kind}{size}"
        if repeat:
            seq = ref[contig]
            anchor = seq[pos] if seq[pos] != repeat[0] else rng.base(exclude=repeat[0])
            ref[contig] = seq[:pos] + anchor + repeat + seq[pos + 1 + len(repeat) :]
            inserted = repeat[:size]
            tag = f"indel_{kind}{size}_in_{repeat[:size]}x{len(repeat) // size}"
        indel = make_variant(rng, ref, contig, pos, kind, size, tag=tag, inserted=inserted)
        variants.append(indel)
        variants.append(
            make_variant(
                rng, ref, contig, snp_pos, "snp", genotypes=list(indel["gts"]), tag="snp_near_indel"
            )
        )
        claim(contig, pos, snp_pos + 1)

    # SNPs that every sample carries on both haplotypes: fail AC=AN.
    for _ in range(3):
        while True:
            contig = ["ctgA", "ctgB"][rng.below(2)]
            pos = rng.between(200, CONTIG_LENGTHS[contig] - 200)
            if in_intervals(contig, pos, 400, blocked) or not free(contig, pos, pos + 1, 30):
                continue
            break
        variants.append(
            make_variant(
                rng,
                ref,
                contig,
                pos,
                "snp",
                genotypes=[(1, 1)] * len(SAMPLES),
                tag="snp_all_homalt",
            )
        )
        claim(contig, pos, pos + 1)

    # SNPs inside the near-identical copies, away from the copy edges (where
    # unique flanking sequence anchors pairs) and from the copy-specific
    # differences: fail MQ<30. The last one is also hom-alt everywhere, so it
    # fails both filters.
    all_homalt = [(1, 1)] * len(SAMPLES)
    near_sites = [
        ("ctgA", NEAR_SRC, 600, None, "snp_low_mq"),
        ("ctgA", NEAR_SRC, 1000, None, "snp_low_mq"),
        ("ctgB", NEAR_DEST, 1400, None, "snp_low_mq"),
        ("ctgB", NEAR_DEST, 800, all_homalt, "snp_all_homalt_low_mq"),
    ]
    for contig, base, offset, genotypes, tag in near_sites:
        offset += rng.between(-20, 20)
        while any(abs(offset - d) < 15 for d in near_diffs):
            offset += 7
        pos = base + offset
        gts = genotypes if genotypes is not None else het_genotypes(rng)
        variants.append(make_variant(rng, ref, contig, pos, "snp", genotypes=gts, tag=tag))
        claim(contig, pos, pos + 1)

    # Ordinary SNPs in unique sequence.
    n_plain = 0
    while n_plain < 34:
        contig = ["ctgA", "ctgA", "ctgB", "ctgC"][rng.below(4)]
        pos = rng.between(200, CONTIG_LENGTHS[contig] - 200)
        if in_intervals(contig, pos, 200, blocked) or not free(contig, pos, pos + 1, 30):
            continue
        variants.append(make_variant(rng, ref, contig, pos, "snp", tag="snp"))
        claim(contig, pos, pos + 1)
        n_plain += 1

    variants.sort(key=lambda v: (list(CONTIG_LENGTHS).index(v["contig"]), v["pos"]))
    return variants


def build_haplotypes(ref, variants):
    """Return {sample: [{contig: seq}, {contig: seq}]}."""
    haplotypes = {}
    for s_idx, sample in enumerate(SAMPLES):
        haps = []
        for h in (0, 1):
            contigs = {}
            for contig, seq in ref.items():
                pieces = []
                cursor = 0
                for v in variants:
                    if v["contig"] != contig or v["gts"][s_idx][h] == 0:
                        continue
                    pieces.append(seq[cursor : v["pos"]])
                    pieces.append(v["alt"])
                    cursor = v["pos"] + len(v["ref"])
                pieces.append(seq[cursor:])
                contigs[contig] = "".join(pieces)
            haps.append(contigs)
        haplotypes[sample] = haps
    return haplotypes


def binned_quality(rng, weights):
    """Draw from NovaSeq-style quality bins, given (quality, weight) pairs."""
    r = rng.random()
    for q, w in weights:
        if r < w:
            return q
        r -= w
    return weights[-1][0]


# Binned qualities keep the fastqs small and match modern Illumina output.
QUAL_BODY = [(37, 0.90), (23, 0.08), (12, 0.02)]
QUAL_TAIL = [(37, 0.75), (23, 0.20), (12, 0.05)]
QUAL_ERROR = [(12, 0.5), (23, 0.3), (2, 0.2)]


def quality_and_errors(rng, seq):
    """Apply substitution errors and return (seq, qual)."""
    out = []
    qual = []
    for i, base in enumerate(seq):
        if rng.random() < ERROR_RATE:
            out.append(rng.base(exclude=base))
            q = binned_quality(rng, QUAL_ERROR)
        else:
            out.append(base)
            q = binned_quality(rng, QUAL_TAIL if i >= READ_LEN - 20 else QUAL_BODY)
        qual.append(chr(q + 33))
    return "".join(out), "".join(qual)


def pick_fragment(rng, haps, length):
    h = rng.below(2)
    contigs = haps[h]
    total = sum(len(s) for s in contigs.values())
    r = rng.below(total)
    for seq in contigs.values():
        if r < len(seq):
            break
        r -= len(seq)
    start = rng.below(len(seq) - length + 1)
    frag = seq[start : start + length]
    if rng.random() < 0.5:
        frag = revcomp(frag)
    return frag


def insert_length(rng):
    while True:
        length = int(round(rng.gauss(INSERT_MEAN, INSERT_SD)))
        if INSERT_MIN <= length <= INSERT_MAX:
            return length


def reads_from_fragment(rng, frag):
    r1 = frag[:READ_LEN]
    r2 = revcomp(frag)[:READ_LEN]
    return r1, r2


def short_insert_reads(rng, frag):
    r1 = frag + ADAPTER_R1
    r2 = revcomp(frag) + ADAPTER_R2
    r1 = (r1 + rng.seq(READ_LEN))[:READ_LEN]
    r2 = (r2 + rng.seq(READ_LEN))[:READ_LEN]
    return r1, r2


def simulate_sample(rng, haps):
    n_dup = int(round(PAIRS_PER_SAMPLE * DUPLICATE_FRACTION))
    n_short = int(round(PAIRS_PER_SAMPLE * SHORT_INSERT_FRACTION))
    n_junk = int(round(PAIRS_PER_SAMPLE * JUNK_MATE_FRACTION))
    n_normal = PAIRS_PER_SAMPLE - n_dup - n_short - n_junk

    fragments = [pick_fragment(rng, haps, insert_length(rng)) for _ in range(n_normal)]
    pairs = [reads_from_fragment(rng, f) for f in fragments]
    # PCR duplicates: the same fragment read again, with fresh errors.
    for _ in range(n_dup):
        pairs.append(reads_from_fragment(rng, fragments[rng.below(len(fragments))]))
    for _ in range(n_short):
        frag = pick_fragment(rng, haps, rng.between(SHORT_INSERT_MIN, SHORT_INSERT_MAX))
        pairs.append(short_insert_reads(rng, frag))
    for _ in range(n_junk):
        r1, _ = reads_from_fragment(rng, pick_fragment(rng, haps, insert_length(rng)))
        pairs.append((r1, rng.seq(READ_LEN)))

    rng.shuffle(pairs)
    out = []
    for r1, r2 in pairs:
        s1, q1 = quality_and_errors(rng, r1)
        s2, q2 = quality_and_errors(rng, r2)
        out.append((s1, q1, s2, q2))
    return out


def write_gzip(path, text):
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(text.encode())


def write_fasta(path, ref):
    lines = []
    for contig, seq in ref.items():
        lines.append(f">{contig}")
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    Path(path).write_text("\n".join(lines) + "\n")


def write_gff(path):
    lines = ["##gff-version 3"]
    for i, (contig, start, end, strand) in enumerate(GENES, 1):
        lines.append(
            f"{contig}\tsim\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene{i};Name=gene{i}"
        )
    Path(path).write_text("\n".join(lines) + "\n")


def write_truth(path, variants):
    header = ["contig", "pos", "ref", "alt", "tag", *SAMPLES]
    lines = ["\t".join(header)]
    for v in variants:
        gts = [f"{a}|{b}" for a, b in v["gts"]]
        # VCF-style 1-based position of the first REF base.
        lines.append(
            "\t".join([v["contig"], str(v["pos"] + 1), v["ref"], v["alt"], v["tag"], *gts])
        )
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory (default: this script's directory)",
    )
    args = parser.parse_args()

    rng = Rng(SEED)
    ref, near_diffs = build_reference(rng)
    variants = place_variants(rng, ref, near_diffs)
    haplotypes = build_haplotypes(ref, variants)

    outdir = args.outdir
    (outdir / "fastq").mkdir(parents=True, exist_ok=True)
    write_fasta(outdir / "reference.fasta", ref)
    write_gff(outdir / "genes.gff")
    write_truth(outdir / "truth.tsv", variants)

    for sample in SAMPLES:
        pairs = simulate_sample(rng, haplotypes[sample])
        r1 = []
        r2 = []
        for n, (s1, q1, s2, q2) in enumerate(pairs, 1):
            r1.append(f"@SIM.{n}\n{s1}\n+\n{q1}\n")
            r2.append(f"@SIM.{n}\n{s2}\n+\n{q2}\n")
        write_gzip(outdir / "fastq" / f"{sample}_1.fastq.gz", "".join(r1))
        write_gzip(outdir / "fastq" / f"{sample}_2.fastq.gz", "".join(r2))


if __name__ == "__main__":
    main()
