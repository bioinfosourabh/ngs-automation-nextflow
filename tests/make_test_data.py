#!/usr/bin/env python3
"""Generate the small synthetic data set used by the `test` profile.

The data are entirely simulated: a random two-contig reference, germline
variants for three patients, somatic SNVs in two tumours, and paired-end
reads simulated with `wgsim` (from samtools). Resource VCFs (dbSNP-like known
sites, an allele-frequency-only germline resource, a panel of normals and a
common-variant set for contamination estimation) are written to match.

Requirements: python >= 3.8, `wgsim`, `bgzip` and `tabix` on PATH.

    python3 tests/make_test_data.py tests/data
"""
import gzip
import os
import random
import shutil
import subprocess
import sys

SEED = 20260101
CONTIGS = {"chrT1": 30000, "chrT2": 10000}
READ_LEN = 100
COVERAGE = 30


def rand_seq(rng, n, gc=0.41):
    at = (1 - gc) / 2
    gcw = gc / 2
    return "".join(rng.choices("ACGT", weights=[at, gcw, gcw, at], k=n))


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def pick_sites(rng, ref, n, exclude, margin=300):
    sites = []
    while len(sites) < n:
        contig = rng.choice(list(ref))
        pos = rng.randint(margin, len(ref[contig]) - margin)  # 1-based
        if any(c == contig and abs(p - pos) < 30 for c, p in exclude + sites):
            continue
        sites.append((contig, pos))
    return sorted(sites)


def make_variant(rng, ref, contig, pos, kind):
    refbase = ref[contig][pos - 1]
    if kind == "snv":
        alt = rng.choice([b for b in "ACGT" if b != refbase])
        return (contig, pos, refbase, alt)
    if kind == "del":
        length = rng.randint(1, 4)
        return (contig, pos, ref[contig][pos - 1:pos + length], refbase)
    ins = rand_seq(rng, rng.randint(1, 4))
    return (contig, pos, refbase, refbase + ins)


def apply_variants(ref, variants):
    """Return a haplotype with the given (non-overlapping) variants applied."""
    hap = {}
    for contig, seq in ref.items():
        out, last = [], 0
        for c, pos, r, a in sorted(v for v in variants if v[0] == contig):
            out.append(seq[last:pos - 1])
            out.append(a)
            last = pos - 1 + len(r)
        out.append(seq[last:])
        hap[contig] = "".join(out)
    return hap


def simulate(hap_fa, n_pairs, seed, prefix, r1, r2):
    tmp1, tmp2 = r1 + ".tmp", r2 + ".tmp"
    subprocess.run(
        ["wgsim", "-N", str(n_pairs), "-1", str(READ_LEN), "-2", str(READ_LEN),
         "-e", "0.001", "-r", "0", "-R", "0", "-d", "350", "-s", "40",
         "-S", str(seed), hap_fa, tmp1, tmp2],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for tmp, final in ((tmp1, r1), (tmp2, r2)):
        with open(tmp) as src, open(final, "a") as dst:
            for i, line in enumerate(src):
                if i % 4 == 0:
                    line = "@" + prefix + "_" + line[1:]
                dst.write(line)
        os.remove(tmp)


def gzip_file(path):
    with open(path, "rb") as src, gzip.GzipFile(path + ".gz", "wb", mtime=0) as dst:
        shutil.copyfileobj(src, dst)
    os.remove(path)


def write_vcf(path, records, info_header="", extra_info=None):
    header = ["##fileformat=VCFv4.2"]
    for name, length in CONTIGS.items():
        header.append(f"##contig=<ID={name},length={length}>")
    if info_header:
        header.append(info_header)
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    with open(path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        for i, (c, p, r, a) in enumerate(sorted(records, key=lambda v: (list(CONTIGS).index(v[0]), v[1]))):
            info = extra_info(i) if extra_info else "."
            fh.write(f"{c}\t{p}\trs{9000000 + i}\t{r}\t{a}\t.\tPASS\t{info}\n")
    subprocess.run(["bgzip", "-f", path], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", path + ".gz"], check=True)


def main(outdir):
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)
    ref = {name: rand_seq(rng, length) for name, length in CONTIGS.items()}
    write_fasta(os.path.join(outdir, "genome.fasta"), ref)

    # Germline variants: a shared common set plus private variants per patient.
    kinds = ["snv"] * 8 + ["del", "ins"]
    used = []
    common_sites = pick_sites(rng, ref, 30, used)
    used += common_sites
    common = [make_variant(rng, ref, c, p, rng.choice(kinds)) for c, p in common_sites]
    patients = {}
    for patient in ("P1", "P2", "P3"):
        private_sites = pick_sites(rng, ref, 8, used)
        used += private_sites
        private = [make_variant(rng, ref, c, p, rng.choice(kinds)) for c, p in private_sites]
        carried = [v for v in common if rng.random() < 0.6] + private
        # Split into two haplotypes: about one third homozygous.
        hap1, hap2 = [], []
        for v in carried:
            r = rng.random()
            if r < 0.33:
                hap1.append(v); hap2.append(v)
            elif r < 0.66:
                hap1.append(v)
            else:
                hap2.append(v)
        patients[patient] = (hap1, hap2)

    somatic_sites = pick_sites(rng, ref, 12, used)
    used += somatic_sites
    somatic = [make_variant(rng, ref, c, p, "snv") for c, p in somatic_sites]

    fastq_dir = os.path.join(outdir, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)
    total_len = sum(CONTIGS.values())
    pairs_total = total_len * COVERAGE // (2 * READ_LEN)

    def make_sample(sample, components, lanes=1):
        """components: list of (haplotype_variants, fraction)."""
        for lane in range(1, lanes + 1):
            r1 = os.path.join(fastq_dir, f"{sample}_L00{lane}_R1.fastq")
            r2 = os.path.join(fastq_dir, f"{sample}_L00{lane}_R2.fastq")
            for f in (r1, r2):
                if os.path.exists(f):
                    os.remove(f)
            for k, (variants, frac) in enumerate(components):
                hap_fa = os.path.join(outdir, "hap.tmp.fa")
                write_fasta(hap_fa, apply_variants(ref, variants))
                n = int(pairs_total * frac / lanes)
                seed = SEED + lane * 100 + k * 7 + sum(map(ord, sample))
                simulate(hap_fa, n, seed, f"{sample}_L{lane}_h{k}", r1, r2)
                os.remove(hap_fa)
            gzip_file(r1)
            gzip_file(r2)

    h1, h2 = patients["P1"]
    make_sample("NORMAL1", [(h1, 0.5), (h2, 0.5)], lanes=2)
    make_sample("TUMOUR1", [(h1, 0.35), (h2, 0.35), (h1 + somatic[:8], 0.30)])
    h1, h2 = patients["P2"]
    make_sample("NORMAL2", [(h1, 0.5), (h2, 0.5)])
    h1, h2 = patients["P3"]
    make_sample("TUMOUR3", [(h1, 0.30), (h2, 0.30), (h2 + somatic[6:], 0.40)])

    germline_all = sorted(set(common + [v for p in patients.values() for h in p for v in h]))
    snvs = [v for v in germline_all if len(v[2]) == 1 and len(v[3]) == 1]
    indels = [v for v in germline_all if v not in snvs]

    af_hdr = '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    write_vcf(os.path.join(outdir, "dbsnp.vcf"), snvs[::2] + indels)
    write_vcf(os.path.join(outdir, "known_indels.vcf"), indels)
    write_vcf(os.path.join(outdir, "germline_resource.vcf"), germline_all, af_hdr,
              lambda i: f"AF={0.05 + (i % 7) * 0.05:.2f}")
    write_vcf(os.path.join(outdir, "common_variants.vcf"), snvs, af_hdr,
              lambda i: f"AF={0.05 + (i % 4) * 0.04:.2f}")
    pon_site = pick_sites(rng, ref, 1, used)[0]
    write_vcf(os.path.join(outdir, "pon.vcf"), [make_variant(rng, ref, *pon_site, "snv")])

    with open(os.path.join(outdir, "targets.bed"), "w") as fh:
        fh.write("chrT1\t150\t14000\nchrT1\t15000\t29850\nchrT2\t150\t9850\n")

    # Truth sets, for checking a real (non-stub) test run by eye.
    with open(os.path.join(outdir, "truth_somatic.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tpatient\n")
        for v in somatic[:8]:
            fh.write("\t".join(map(str, v)) + "\tP1\n")
        for v in somatic[6:]:
            fh.write("\t".join(map(str, v)) + "\tP3\n")
    with open(os.path.join(outdir, "truth_germline.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tpatient\tgenotype\n")
        for patient, (hap1, hap2) in patients.items():
            for v in sorted(set(hap1 + hap2)):
                gt = "1/1" if (v in hap1 and v in hap2) else "0/1"
                fh.write("\t".join(map(str, v)) + f"\t{patient}\t{gt}\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
