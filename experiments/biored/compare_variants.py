"""Paired comparison of two single-shot variants on the same BioRED abstracts (e.g. the card vs the minimal prompt).

Usage: python experiments/biored/compare_variants.py [variant_a] [variant_b] [run_set]
  defaults: v1_card_lean_medium v6_minimal_lean_medium biored_all

Scores only abstracts both variants finished, per BioRED split and pooled, at three operating points:
strict (supported), routed (supported + uncertain whose reasons are only hedged_claim and/or scope) and
supported_or_uncertain. The difference B - A gets a paired bootstrap 95% interval (abstracts resampled once for
both variants, 2,000 draws, fixed seed) and the share of draws in which B's F1 is not lower.
Also reports per-variant cost, first-try validity and prediction counts. Writes runs/<run_set>/<b>/compare_<a>.json.
"""
import csv
import json
import random
import sys

import evaluate
from evaluate import HERE, RUNS, prf, score_paper, total  # prf returns (precision, recall, f1)

ROUTED = {"hedged_claim", "scope"}
_read = evaluate.read


def read_routed(path):
    """pairs.csv with uncertain pairs whose reasons are only hedged_claim/scope promoted to supported."""
    rows = _read(path)
    if path.name == "pairs.csv":
        for r in rows:
            reasons = set(json.loads(r["unresolved_reason"] or "[]"))
            if r["status"] == "uncertain" and reasons and reasons <= ROUTED:
                r["status"] = "supported"
    return rows


def scores(variant, run_set, papers, gold):
    base = RUNS / run_set / variant
    out = {"strict": {}, "routed": {}, "supported_or_uncertain": {}}
    for p in papers:
        out["strict"][p] = score_paper(p, gold[p], base / p / "outputs")[0]
        out["supported_or_uncertain"][p] = score_paper(p, gold[p], base / p / "outputs", ("supported", "uncertain"))[0]
        evaluate.read = read_routed
        try:
            out["routed"][p] = score_paper(p, gold[p], base / p / "outputs")[0]
        finally:
            evaluate.read = _read
    return out


def run_stats(variant, run_set, papers):
    with (RUNS / run_set / variant / "summary.csv").open(encoding="utf-8-sig") as stream:
        rows = {r["paper_id"]: r for r in csv.DictReader(stream)}
    rows = [rows[p] for p in papers if p in rows]
    return {"cost_usd": round(sum(float(r["cost_usd"]) for r in rows), 2),
            "first_try": f"{sum(r['calls'] == '1' for r in rows)}/{len(rows)}",
            "supported": sum(int(r["supported"] or 0) for r in rows),
            "uncertain": sum(int(r["uncertain"] or 0) for r in rows),
            "genes": sum(int(r["genes"] or 0) for r in rows), "targets": sum(int(r["targets"] or 0) for r in rows)}


def prf_dict(counts):
    return dict(zip(("precision", "recall", "f1"), prf(total(counts))))


def paired(a, b, papers, rng):
    pa, pb = prf_dict([a[p] for p in papers]), prf_dict([b[p] for p in papers])
    diffs = {k: [] for k in ("precision", "recall", "f1")}
    for _ in range(2000):
        sample = [rng.choice(papers) for _ in papers]
        sa, sb = prf_dict([a[p] for p in sample]), prf_dict([b[p] for p in sample])
        for k in diffs:
            diffs[k].append(sb[k] - sa[k])
    result = {"a": pa, "b": pb}
    for k, d in diffs.items():
        d.sort()
        result[f"diff_{k}"] = {"value": round(pb[k] - pa[k], 3), "ci95": [round(d[49], 3), round(d[1949], 3)]}
    result["share_b_f1_not_lower"] = round(sum(x >= 0 for x in diffs["f1"]) / len(diffs["f1"]), 3)
    return result


def main(a="v1_card_lean_medium", b="v6_minimal_lean_medium", run_set="biored_all"):
    gold = json.loads((HERE / "gold_all.json").read_text(encoding="utf-8"))
    done = lambda v, p: (RUNS / run_set / v / p / "outputs" / "validation.json").exists()
    papers = sorted(p for p in gold if done(a, p) and done(b, p))
    split = {p["paper_id"]: p["split"] for p in json.loads((HERE / "sample_all.json").read_text(encoding="utf-8"))["papers"]}
    sa, sb = scores(a, run_set, papers, gold), scores(b, run_set, papers, gold)
    rng = random.Random(20261028)
    groups = {name: [p for p in papers if split[p] == name] for name in ("Train", "Dev", "Test")}
    groups["all scored"] = papers
    result = {"a": a, "b": b, "run_set": run_set, "papers": len(papers),
              "runs": {a: run_stats(a, run_set, papers), b: run_stats(b, run_set, papers)}, "groups": {}}
    for name, subset in groups.items():
        if subset:
            result["groups"][name] = {"papers": len(subset),
                                      **{point: paired(sa[point], sb[point], subset, rng) for point in sa}}
    path = RUNS / run_set / b / f"compare_{a}.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{a} (A) vs {b} (B) on {len(papers)} abstracts; written to {path}")
    for v, s in result["runs"].items():
        print(f"  {v}: " + ", ".join(f"{k} {x}" for k, x in s.items()))
    for name, g in result["groups"].items():
        print(f"\n{name} ({g['papers']} abstracts)")
        for point in sa:
            r = g[point]
            print(f"  {point:23s} A P/R/F1 {r['a']['precision']:.3f}/{r['a']['recall']:.3f}/{r['a']['f1']:.3f}"
                  f"  B {r['b']['precision']:.3f}/{r['b']['recall']:.3f}/{r['b']['f1']:.3f}"
                  f"  dF1 {r['diff_f1']['value']:+.3f} {r['diff_f1']['ci95']}"
                  f"  dP {r['diff_precision']['value']:+.3f} {r['diff_precision']['ci95']}"
                  f"  dR {r['diff_recall']['value']:+.3f} {r['diff_recall']['ci95']}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
