#!/usr/bin/env python3
"""
ATB cross-species FILTER -- the pass/fail half of the pipeline-workflow replacement for
the bg_excl set-diff step (MARKER_VALIDATION_PLAN.md Tier B). Reads the JSONL produced by
run_atb_pseudoalign.sh (themisto2 threshold-pseudoalign --threshold 0.01 --denominator all
--report-hit-counts) and, per marker:

    hit_frac[colour] = kmer_hits[colour] / marker_kmer_count      (marker_kmer_count = len - k + 1)

    PASS  iff  hit_frac[target_species] >= --min-within   (default 0.95)
          AND  every non-target, non-excluded colour's hit_frac <= --max-outside (default 0.05)
    FLAG  otherwise (hits the target but also leaks into >=1 off-target species above tolerance)
    ABSENT if the target species is never hit at all (marker doesn't come from this species'
           index in the first place -- almost certainly a Tier A problem, not a Tier B one)

This is the SAME >=95%/<=5% tolerance lineage_specificity_filter.py's lineage_view() and
Tier A already use, just applied to ATB species colours instead of within-species lineage
colours -- see the conversation note: --min-within is enforced here (unlike the current
prototype validate.py tier-b, which treats it as vestigial and only gates on tf>0). If you
want to match validate.py's current behaviour instead, pass --min-within 0.

Outputs (all under --out prefix):
    <out>_validation.tsv         one row per marker: fractions, verdict, reason
    <out>_species_detail.tsv     one row per marker x non-target-species hit (uncapped)
    <out>_summary.txt            counts, pass rate, top offending species
    <out>_PASS.fasta             ONLY the markers that PASSed -- this is the actual set-diff
                                  replacement's output: feed this forward instead of the raw
                                  candidate FASTA
    <out>_FLAG.fasta             markers dropped for off-target leakage (kept for inspection,
                                  not forwarded)
    <out>_PASS.ids / _FLAG.ids / _ABSENT.ids
"""
import argparse
import collections
import json
import sys
from pathlib import Path


def read_fasta(path):
    """id -> (header_line_without_gt, sequence) preserving input order."""
    seqs = collections.OrderedDict()
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]  # ID = first whitespace-delimited token
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs


def _clean_atb_name(raw):
    """Strip ATB's 'per_species_unitigs/<name>-unitigs-k31.fna' wrapper down to the bare
    species/colour name. Leaves anything that doesn't match the wrapper untouched (so a
    plain name-per-line file, or a genuinely different naming convention, still works)."""
    name = raw
    prefix = "per_species_unitigs/"
    if name.startswith(prefix):
        name = name[len(prefix):]
    suffix = "-unitigs-k31.fna"
    if name.endswith(suffix):
        name = name[:-len(suffix)]
    return name


def load_color_names(path):
    """color_names.txt -> {colour_id(int): name}. Tolerates 'id\\tname', ATB's actual
    'id\\tper_species_unitigs/<name>-unitigs-k31.fna' format (unwrapped via _clean_atb_name),
    or a plain name-per-line file where the line number (0-based) is the colour id, matching
    themisto2's export convention."""
    names = {}
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    looks_indexed = any("\t" in l for l in lines[:5])
    if looks_indexed:
        for l in lines:
            if not l:
                continue
            cid, name = l.split("\t", 1)
            names[int(cid)] = _clean_atb_name(name)
    else:
        for i, name in enumerate(lines):
            if name:
                names[i] = _clean_atb_name(name)
    return names


def expand_target_species(base_names, names):
    """Resolve user-given base species names (e.g. 'streptococcus_pneumoniae') against the
    full colour-name set, picking up ATB's lettered-split chunks (streptococcus_pneumoniaea,
    ...b, ...c, ...) alongside the base colour itself -- a species is split iff a colour
    equals base + exactly one trailing lowercase letter.

    Returns (resolved, ambiguous):
      resolved  -- {colour_id: (base_name, colour_name)} for the base name itself and any
                   clean lettered chunk of it.
      ambiguous -- [(colour_id, colour_name, base_name)] for colours that merely CONTAIN a
                   base name as a substring without fitting that clean pattern (ATB has at
                   least one glued-together name, e.g.
                   'legionella_pneumophilastreptococcus_pneumoniae' -- two full genus_species
                   names concatenated with no separator, not a real split of either). These
                   can't be confidently called target or non-target, so they're excluded from
                   both sides of the check and reported as a warning instead of silently
                   scored either way.
    """
    resolved = {}
    ambiguous = []
    for base in base_names:
        for cid, name in names.items():
            if name == base:
                resolved[cid] = (base, name)
            elif name.startswith(base) and len(name) == len(base) + 1 and name[-1].isalpha() and name[-1].islower():
                resolved[cid] = (base, name)  # clean lettered split, e.g. ...pneumoniaea
            elif base in name:
                ambiguous.append((cid, name, base))
    return resolved, ambiguous


def parse_jsonl(path):
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                yield json.loads(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jsonl", required=True, help="themisto2 threshold-pseudoalign JSONL")
    ap.add_argument("--fasta", required=True, help="the same marker FASTA used as the query")
    ap.add_argument("--color-names", required=True)
    ap.add_argument("--target-species", required=True, nargs="+",
                     help="one or more ATB colour names counted as on-target (ATB splits big "
                          "species into lettered chunks, e.g. streptococcus_pneumoniae{,a,b,c})")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--min-within", type=float, default=0.95,
                     help="target species hit_frac must be >= this to PASS (default 0.95)")
    ap.add_argument("--max-outside", type=float, default=0.05,
                     help="every non-target, non-excluded species hit_frac must be <= this to "
                          "PASS (default 0.05 -- same tolerance the pipeline's own specificity "
                          "filter uses)")
    ap.add_argument("--exclude-species", nargs="*", default=["unknown"],
                     help="ATB colours dropped from the non-target check entirely (default: "
                          "'unknown', ATB's catch-all bucket)")
    ap.add_argument("--out", required=True, help="output path prefix")
    args = ap.parse_args()

    names = load_color_names(args.color_names)
    seqs = read_fasta(args.fasta)
    lengths = {mid: len(seq) for mid, seq in seqs.items()}
    k = args.k
    tgt_names = set(args.target_species)
    resolved, ambiguous = expand_target_species(tgt_names, names)
    tgt_ids = set(resolved)
    if not tgt_ids:
        sys.exit(f"no target species {sorted(tgt_names)} found in {args.color_names}")
    matched_bases = {base for base, _ in resolved.values()}
    missing_bases = tgt_names - matched_bases
    if missing_bases:
        sys.exit(f"target species not found in {args.color_names}: {sorted(missing_bases)}")
    extra_chunks = sorted(name for base, name in resolved.values() if name not in tgt_names)
    if extra_chunks:
        print(f"[expand_target_species] also treating as target (ATB lettered split): {extra_chunks}",
              file=sys.stderr)
    # Ambiguous colours (contain a target name as a substring but aren't a clean lettered
    # split of it) are excluded from BOTH the target and non-target sides of the check --
    # never guessed either way, only reported.
    ambig_ids = {cid for cid, name, base in ambiguous}
    if ambiguous:
        print(f"[expand_target_species] {len(ambiguous)} colour(s) ambiguously related to a "
              f"target name -- excluded from scoring, see {args.out}_warnings.tsv", file=sys.stderr)
    excl = set(args.exclude_species)  # by name (--exclude-species)

    rows = []
    species_detail = []
    seen = set()
    for rec in parse_jsonl(args.jsonl):
        mid = str(rec["name"])
        seen.add(mid)
        nk = max(lengths.get(mid, 0) - k + 1, 0)
        cols = rec.get("colors", [])
        hits = rec.get("kmer_hits", [])
        if len(hits) != len(cols):
            rows.append(dict(marker_id=mid, n_kmers=nk, target_frac="", max_nontarget_frac="",
                              max_nontarget_species="", verdict="FLAG",
                              reason="malformed record (colors/kmer_hits length mismatch)"))
            continue
        tf = max((h / nk for c, h in zip(cols, hits) if c in tgt_ids), default=0.0) if nk else 0.0

        nt = []  # (frac, hits, name) for every non-target, non-excluded colour hit
        for c, h in zip(cols, hits):
            name = names.get(c, str(c))
            if c in tgt_ids or c in ambig_ids:
                continue
            frac = (h / nk) if nk else 0.0
            species_detail.append(dict(marker_id=mid, species=name, frac=round(frac * 100, 2),
                                        n_hits=h, n_kmers=nk))
            if name in excl:
                continue
            nt.append((frac, h, name))
        nt.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
        mnt, mnts = (nt[0][0], nt[0][2]) if nt else (0.0, "-")

        if tf <= 0:
            verdict, reason = "ABSENT", "target species never hit"
        elif tf < args.min_within:
            verdict, reason = "FLAG", f"target_frac {tf*100:.1f}% < min-within {args.min_within*100:g}%"
        elif mnt > args.max_outside:
            verdict = "FLAG"
            reason = f"{mnts} {mnt*100:.1f}% > max-outside {args.max_outside*100:g}%"
        else:
            verdict, reason = "PASS", "-"

        nt_str = "; ".join(f"{name} {frac*100:.2f}% ({h}/{nk})" for frac, h, name in nt[:20])
        rows.append(dict(marker_id=mid, n_kmers=nk, n_nontarget_species=len(nt),
                          target_frac=round(tf * 100, 2), max_nontarget_frac=round(mnt * 100, 2),
                          max_nontarget_species=mnts, nontarget_species=nt_str,
                          verdict=verdict, reason=reason))

    # markers in the FASTA that themisto never emitted a record for at all (e.g. zero hits
    # anywhere, including target) -- still ABSENT, not silently dropped from the accounting
    for mid in seqs:
        if mid not in seen:
            rows.append(dict(marker_id=mid, n_kmers=max(lengths[mid] - k + 1, 0),
                              n_nontarget_species=0, target_frac=0.0, max_nontarget_frac=0.0,
                              max_nontarget_species="-", nontarget_species="",
                              verdict="ABSENT", reason="no pseudoalign record (zero hits)"))

    out = Path(args.out)
    cols = ["marker_id", "n_kmers", "n_nontarget_species", "target_frac", "max_nontarget_frac",
            "max_nontarget_species", "nontarget_species", "verdict", "reason"]
    with open(out.with_name(out.name + "_validation.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # colours that couldn't be confidently resolved as target or non-target (see
    # expand_target_species) -- excluded from every marker's scoring above; recorded here
    # once, per colour, not per marker, so this is a colour-level warning, not a verdict.
    if ambiguous:
        with open(out.with_name(out.name + "_warnings.tsv"), "w") as fh:
            fh.write("colour_id\tcolour_name\tmatched_target_substring\treason\n")
            for cid, name, base in ambiguous:
                fh.write(f"{cid}\t{name}\t{base}\tcontains target name as a substring but is "
                         f"not a clean lettered split of it -- excluded from scoring\n")

    if species_detail:
        dcols = ["marker_id", "species", "frac", "n_hits", "n_kmers"]
        with open(out.with_name(out.name + "_species_detail.tsv"), "w") as fh:
            fh.write("\t".join(dcols) + "\n")
            for d in species_detail:
                fh.write("\t".join(str(d[c]) for c in dcols) + "\n")

    # --- drop the FASTAs by verdict: this is the actual set-diff replacement's output ---
    by_verdict = collections.defaultdict(list)
    for r in rows:
        by_verdict[r["verdict"]].append(r["marker_id"])
    for verdict in ("PASS", "FLAG", "ABSENT"):
        ids = by_verdict.get(verdict, [])
        fa_path = out.with_name(out.name + f"_{verdict}.fasta")
        with open(fa_path, "w") as fh:
            for mid in ids:
                seq = seqs.get(mid)
                if seq is not None:
                    fh.write(f">{mid}\n{seq}\n")
        out.with_name(out.name + f"_{verdict}.ids").write_text(
            "\n".join(ids) + ("\n" if ids else ""))

    # --- summary report ---
    n = len(rows)
    npass = len(by_verdict.get("PASS", []))
    nflag = len(by_verdict.get("FLAG", []))
    nabsent = len(by_verdict.get("ABSENT", []))
    nscored = n - nabsent
    sc_all = collections.Counter()
    for d in species_detail:
        if d["frac"] > args.max_outside * 100:
            sc_all[d["species"]] += 1
    lines = [
        "tier                 : ATB cross-species check (set-diff replacement)",
        f"target species       : {' '.join(sorted(tgt_names))}"
        + (f"  (+ {len(extra_chunks)} ATB lettered split(s): {extra_chunks})" if extra_chunks else ""),
        f"excluded colours     : {','.join(sorted(excl)) or '-'}"
        + (f"  + {len(ambiguous)} ambiguous colour(s) (see {out.name}_warnings.tsv)" if ambiguous else ""),
        f"verdict rule         : PASS iff target_frac >= {args.min_within*100:g}% "
        f"AND every other species' hit_frac <= {args.max_outside*100:g}%",
        f"markers evaluated    : {n}  ({nabsent} ABSENT, {nscored} scored)",
        f"PASS                 : {npass}" + (f"  ({100*npass/nscored:.1f}% of scored)" if nscored else ""),
        f"FLAG                 : {nflag}" + (f"  ({100*nflag/nscored:.1f}% of scored)" if nscored else ""),
        f"ABSENT               : {nabsent}  (target species never hit -- check Tier A / index choice)",
    ]
    if sc_all:
        lines.append(f"off-target species causing FLAGs, by # markers over {args.max_outside*100:g}% "
                      f"(top 20 of {len(sc_all)}):")
        for s, c in sc_all.most_common(20):
            lines.append(f"      {s:<34} {c}")
    lines.append("")
    lines.append(f"FASTA outputs        : {out.name}_PASS.fasta ({npass} seqs) -- forward this; "
                 f"{out.name}_FLAG.fasta ({nflag}) and {out.name}_ABSENT.fasta ({nabsent}) kept for inspection")
    txt = "\n".join(lines) + "\n"
    out.with_name(out.name + "_summary.txt").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
