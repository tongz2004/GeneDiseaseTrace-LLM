"""Paired comparison of two BioTriplex runs on the articles both finished, with a paper-level bootstrap.

Usage: python experiments/biotriplex/compare.py A B [--run-set biotriplex] [--deployment]
  A and B are [<run_set>/]<variant> (full text) or [<run_set>/]<variant>:abstract; the run set defaults to
  --run-set, e.g.
  python experiments/biotriplex/compare.py v6_minimal_lean_medium pubtator3
  python experiments/biotriplex/compare.py v6_minimal_lean_medium v6_minimal_lean_medium:abstract
  python experiments/biotriplex/compare.py biotriplex_gpt/v6_minimal_lean_medium v6_minimal_lean_medium

Articles are paired by PMCID (BTX and BTXA of one article share it). Reports P/R/F1 of each run, recall on the
relations BioTriplex links only outside the abstract, and A - B with 95% intervals (2,000 draws over articles, seed
20261028). Also counts A's supported pairs that are not gold relations by kind. --deployment compares all articles and
counts an article without a valid result as an empty answer. Writes
runs/<run_set>/compare_<A>_vs_<B>[_deployment].json. Run set biotriplex2 uses the second batch's sample and gold.
"""
import argparse
import collections
import importlib.util
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent / "single_shot" / "runs"
_spec = importlib.util.spec_from_file_location("biored_evaluate", HERE.parent / "biored" / "evaluate.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


EMPTY = {"tp_gold": 0, "fn": 0, "fn_symptom": 0, "our_supported": 0, "tp_ours": 0, "uncertain_or_negative": 0,
         "uncertain_or_negative_in_gold": 0, "tp_ours_many_to_one": 0, "tp_gold_many_to_one": 0}


def counts(run_set, spec, articles, gold, deployment=False):
    spec, _, kind = spec.partition(":")
    run_set, _, variant = spec.rpartition("/") if "/" in spec else (run_set, "", spec)
    out = {}
    for a in articles:
        paper = ("BTXA" if kind == "abstract" else "BTX") + a
        o = RUNS / run_set / variant / paper / "outputs"
        if not ((o / "validation.json").exists() or (o / "source.json").exists()):
            if deployment:
                rel = gold[paper]["relations"]
                out[a] = {**EMPTY, "fn": len(rel), "body_gold": sum(r["body_only"] for r in rel), "body_tp": 0,
                          "fp_kinds": collections.Counter()}
            continue
        c, rows = B.score_paper(paper, gold[paper], o)
        rel = gold[paper]["relations"]
        c["body_gold"] = sum(r["body_only"] for r in rel)
        c["body_tp"] = sum(r["body_only"] and row["outcome"] == "TP" for r, row in zip(rel, rows))
        c["fp_kinds"] = collections.Counter(row["detail"] for row in rows if row["outcome"] == "FP")
        out[a] = c
    return out


def metrics(cs):
    t = B.total([{k: v for k, v in c.items() if k != "fp_kinds"} for c in cs])
    p, r, f = B.prf(t)
    body = sum(c["body_tp"] for c in cs) / max(1, sum(c["body_gold"] for c in cs))
    return {"precision": p, "recall": r, "f1": f, "recall_outside_abstract": body}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--run-set", default="biotriplex")
    parser.add_argument("--deployment", action="store_true")
    args = parser.parse_args()
    suffix = "2" if args.run_set == "biotriplex2" else ""
    gold = json.loads((HERE / f"gold_biotriplex{suffix}.json").read_text(encoding="utf-8"))
    articles = [p.removeprefix("BTX") for p in json.loads((HERE / f"sample_biotriplex{suffix}.json").read_text(encoding="utf-8"))["evaluation"]]
    ca = counts(args.run_set, args.a, articles, gold, args.deployment)
    cb = counts(args.run_set, args.b, articles, gold, args.deployment)
    common = sorted(set(ca) & set(cb))
    ma, mb = metrics([ca[a] for a in common]), metrics([cb[a] for a in common])
    rng, diffs = random.Random(20261028), collections.defaultdict(list)
    for _ in range(2000):
        draw = [rng.choice(common) for _ in common]
        da, db = metrics([ca[a] for a in draw]), metrics([cb[a] for a in draw])
        for k in da:
            diffs[k].append(da[k] - db[k])
    ci = {k: [round(sorted(v)[49], 3), round(sorted(v)[1949], 3)] for k, v in diffs.items()}
    result = {"a": args.a, "b": args.b, "articles_compared": len(common),
              "not_in_both": sorted(set(articles) - set(common)),
              "a_metrics": {k: round(v, 3) for k, v in ma.items()}, "b_metrics": {k: round(v, 3) for k, v in mb.items()},
              "a_minus_b": {k: {"diff": round(ma[k] - mb[k], 3), "ci95": ci[k]} for k in ma},
              "a_supported_not_gold": dict(sum((ca[a]["fp_kinds"] for a in common), collections.Counter()))}
    name = f"compare_{args.a}_vs_{args.b}{'_deployment' if args.deployment else ''}".replace(":", "-").replace("/", "-")
    (RUNS / args.run_set / f"{name}.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
