#!/usr/bin/env python3
"""
Collect per-sample QC metrics and apply pass/fail gates.

Inputs (any may be missing → metric reported as NA, gate skipped):
  --mosdepth   <sample>.mosdepth.summary.txt   (+ optional thresholds.bed.gz for WES)
  --markdup    <sample>.markdup.txt            (samtools markdup -f)
  --selfsm     <sample>.selfSM                 (VerifyBamID2)
  --somalier   somalier.samples.tsv            (somalier relate)
Outputs:
  <sample>.qc_summary.tsv       one row, tab-separated, with a final STATUS column
  <sample>_qc_mqc.tsv           MultiQC custom-content table
"""
import argparse, csv, gzip, sys


def read_mosdepth_summary(path):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["chrom"]] = row
    return out


def read_markdup(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if ":" in line:
                k, v = line.rstrip("\n").split(":", 1)
                try:
                    d[k.strip()] = float(v.strip())
                except ValueError:
                    pass
    read = d.get("READ", 0) or 0
    dup = d.get("DUPLICATE TOTAL", 0) or 0
    return (dup / read) if read else None


def read_selfsm(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return float(rows[0]["FREEMIX"]) if rows else None


def read_somalier(path, sample):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("sample_id") == sample:
                # original_pedigree_sex is what we passed in; sex is inferred (1=male, 2=female, 0=unknown)
                inferred = {"1": "male", "2": "female"}.get(row.get("sex", "0"), "unknown")
                return inferred, float(row.get("X_het", "nan") or "nan"), float(row.get("Y_depth_mean", "nan") or "nan")
    return None, None, None


def wes_target_fraction(thresholds_path, min_cov):
    """Fraction of target bases ≥ min_cov from mosdepth --thresholds (needs a column named str(min_cov))."""
    total = covered = 0
    opener = gzip.open if thresholds_path.endswith(".gz") else open
    with opener(thresholds_path, "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            col = header.index(f"{min_cov}X")
        except ValueError:
            return None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            length = int(f[2]) - int(f[1])
            total += length
            covered += int(f[col])
    return covered / total if total else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--seq-type", required=True, choices=["wgs", "wes"])
    ap.add_argument("--expected-sex", default="unknown")
    ap.add_argument("--mosdepth")
    ap.add_argument("--thresholds")
    ap.add_argument("--markdup")
    ap.add_argument("--selfsm")
    ap.add_argument("--somalier")
    ap.add_argument("--max-contamination", type=float, default=0.03)
    ap.add_argument("--min-mean-cov-wgs", type=float, default=25)
    ap.add_argument("--min-target-cov-wes", type=int, default=20)
    ap.add_argument("--min-target-frac-wes", type=float, default=0.90)
    ap.add_argument("--max-duplication", type=float, default=0.25)
    a = ap.parse_args()

    m = {"sample": a.sample, "seq_type": a.seq_type, "expected_sex": a.expected_sex}
    fails, warns = [], []

    # Coverage
    if a.mosdepth:
        s = read_mosdepth_summary(a.mosdepth)
        key = "total_region" if "total_region" in s else "total"
        m["mean_coverage"] = float(s[key]["mean"])
        if a.seq_type == "wgs" and m["mean_coverage"] < a.min_mean_cov_wgs:
            fails.append(f"mean_coverage<{a.min_mean_cov_wgs}")
        if "chrM" in s:
            m["chrM_mean_coverage"] = float(s["chrM"]["mean"])
    if a.seq_type == "wes" and a.thresholds:
        frac = wes_target_fraction(a.thresholds, a.min_target_cov_wes)
        m[f"target_frac_ge{a.min_target_cov_wes}x"] = frac
        if frac is not None and frac < a.min_target_frac_wes:
            fails.append(f"target_frac_ge{a.min_target_cov_wes}x<{a.min_target_frac_wes}")

    # Duplication
    if a.markdup:
        dup = read_markdup(a.markdup)
        m["duplication_rate"] = dup
        if dup is not None and dup > a.max_duplication:
            warns.append(f"duplication>{a.max_duplication}")

    # Contamination
    if a.selfsm:
        fm = read_selfsm(a.selfsm)
        m["freemix"] = fm
        if fm is not None and fm > a.max_contamination:
            fails.append(f"contamination>{a.max_contamination}")

    # Sex
    if a.somalier:
        inferred, xhet, ydepth = read_somalier(a.somalier, a.sample)
        m["inferred_sex"] = inferred
        m["X_het"] = xhet
        m["Y_depth_mean"] = ydepth
        if inferred and a.expected_sex != "unknown" and inferred != "unknown" and inferred != a.expected_sex:
            fails.append(f"sex_mismatch(expected={a.expected_sex},inferred={inferred})")

    m["warnings"] = ";".join(warns) or "."
    m["failures"] = ";".join(fails) or "."
    m["STATUS"] = "FAIL" if fails else ("WARN" if warns else "PASS")

    cols = list(m.keys())
    with open(f"{a.sample}.qc_summary.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join("NA" if m[c] is None else (f"{m[c]:.4f}" if isinstance(m[c], float) else str(m[c])) for c in cols) + "\n")

    # MultiQC custom content
    with open(f"{a.sample}_qc_mqc.tsv", "w") as fh:
        fh.write("# id: fatqvcf_qc_gates\n# section_name: 'fatq-vcf QC gates'\n# plot_type: 'table'\n")
        fh.write("# description: 'Sample-level pass/fail gates (coverage, duplication, contamination, sex). FAIL = do not report.'\n")
        fh.write("Sample\t" + "\t".join(c for c in cols if c != "sample") + "\n")
        fh.write(a.sample + "\t" + "\t".join("NA" if m[c] is None else (f"{m[c]:.4f}" if isinstance(m[c], float) else str(m[c])) for c in cols if c != "sample") + "\n")

    print(f"[qc_gates] {a.sample}: {m['STATUS']}  failures={m['failures']}  warnings={m['warnings']}", file=sys.stderr)


if __name__ == "__main__":
    main()
