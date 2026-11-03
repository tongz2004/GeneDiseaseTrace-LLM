"""Preliminary reproducibility check: the same variant run twice on the same BioRED abstracts.

Usage: python experiments/biored/reproducibility.py [variant ...] [--first biored_all] [--second biored_all_rep2]
  default variants: v1_card_lean_medium v6_minimal_lean_medium

Two runs show whether results move noticeably between calls; they do not estimate run-to-run variance.
Per variant, on the abstracts both runs finished:
- each run's P/R/F1 (one-to-one scorer) and its failures;
- pairs are aligned across runs by mention overlap of both ends (entity and pair IDs are not stable between runs),
  one-to-one; a pair one run lists and the other does not is 'absent' there (unlisted pairs are not_found);
- Jaccard of the two supported sets, status transition counts, and agreement on pairs with a relation status
  (supported / uncertain / explicit_negative) in at least one run. The whole gene x disease matrix is not used,
  because the many not_found pairs would inflate agreement.
Writes runs/<second>/reproducibility.json.
"""
import argparse
import json
from collections import Counter

from adjudication import same, system_pairs
from evaluate import HERE, RUNS, one_to_one, prf, score_paper, total

RELATION = ("supported", "uncertain", "explicit_negative")


def align(first, second):
    """[(status in first, status in second)] over pairs listed with a relation status in either run."""
    matched = one_to_one([[i for i, q in enumerate(second) if same(p, q)] for p in first])
    out = [(first[i]["status"], second[j]["status"]) for i, j in matched.items()]
    out += [(p["status"], "absent") for i, p in enumerate(first) if i not in matched]
    out += [("absent", q["status"]) for j, q in enumerate(second) if j not in set(matched.values())]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("variants", nargs="*", default=["v1_card_lean_medium", "v6_minimal_lean_medium"])
    parser.add_argument("--first", default="biored_all")
    parser.add_argument("--second", default="biored_all_rep2")
    args = parser.parse_args()
    gold = json.loads((HERE / "gold_all.json").read_text(encoding="utf-8"))
    report = {}
    for v in args.variants:
        runs = [RUNS / args.first / v, RUNS / args.second / v]
        attempted = sorted(p.name for p in runs[1].iterdir() if p.name in gold and (p / "run.json").exists())
        ok = lambda base, p: (base / p / "outputs" / "validation.json").exists()
        papers = [p for p in attempted if ok(runs[0], p) and ok(runs[1], p)]
        per_run = []
        for base in runs:
            counts = [score_paper(p, gold[p], base / p / "outputs")[0] for p in papers]
            per_run.append(dict(zip(("precision", "recall", "f1"), (round(x, 3) for x in prf(total(counts))))))
        transitions, related = Counter(), []
        for p in papers:
            pairs = [system_pairs(base / p / "outputs", RELATION) for base in runs]
            for a, b in align(*pairs):
                transitions[f"{a} -> {b}"] += 1
                related.append(a == b)
        sup_both = transitions["supported -> supported"]
        sup_either = sum(n for k, n in transitions.items() if "supported" in k.split(" -> "))
        report[v] = {"papers_compared": len(papers),
                     "failed_or_missing": {args.first: [p for p in attempted if not ok(runs[0], p)],
                                           args.second: [p for p in attempted if not ok(runs[1], p)]},
                     "per_run": dict(zip((args.first, args.second), per_run)),
                     "supported_jaccard": round(sup_both / sup_either, 3) if sup_either else None,
                     "agreement_on_relation_pairs": round(sum(related) / len(related), 3) if related else None,
                     "relation_pairs": len(related), "transitions": dict(transitions.most_common())}
        r = report[v]
        print(f"{v}: {len(papers)} abstracts; P/R/F1 {per_run[0]} vs {per_run[1]}; supported Jaccard {r['supported_jaccard']}; "
              f"agreement on {len(related)} relation pairs {r['agreement_on_relation_pairs']}")
        print("  " + "; ".join(f"{k} {n}" for k, n in transitions.most_common(8)))
    (RUNS / args.second / "reproducibility.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
