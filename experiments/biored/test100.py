"""BioRED Test as a whole: the 73 abstracts with a gene-disease relation (run set biored_all) plus the 27 without
(run set biored_testrest), scored one-to-one with bootstrap 95% intervals over abstracts.

Usage: python experiments/biored/test100.py [variant ...]    (default v1_card_lean_medium v6_minimal_lean_medium pubtator3)

Two accounts per variant:
- valid only: abstracts whose run produced a valid result;
- deployment: every attempted abstract, a failed or missing run counting as an empty answer (its gold relations
  missed, no predictions). A supported pair on the 27 abstracts is a false positive relative to the gold standard.
Also paired differences between every two variants on the abstracts both have.
Writes runs/biored_testrest/test100.json.
"""
import json
import random
import sys

from evaluate import HERE, RUNS, prf, score_paper, total, with_ci

EMPTY = {"tp_gold": 0, "fn": 0, "fn_symptom": 0, "our_supported": 0, "tp_ours": 0, "uncertain_or_negative": 0,
         "uncertain_or_negative_in_gold": 0, "tp_ours_many_to_one": 0, "tp_gold_many_to_one": 0}


def paired(a, b, counts, failed):
    """a - b on the abstracts both have, 2,000 bootstrap draws over abstracts (seed 20261028)."""
    common = sorted(set(counts[a]) - set(failed[a]) - set(failed[b]))
    rng, diffs = random.Random(20261028), []
    point = [x - y for x, y in zip(prf(total([counts[a][p] for p in common])), prf(total([counts[b][p] for p in common])))]
    for _ in range(2000):
        draw = [rng.choice(common) for _ in common]
        diffs.append([x - y for x, y in zip(prf(total([counts[a][p] for p in draw])), prf(total([counts[b][p] for p in draw])))])
    ci = lambda i: [round(sorted(d[i] for d in diffs)[49], 3), round(sorted(d[i] for d in diffs)[1949], 3)]
    return {"abstracts": len(common), **{k: {"diff": round(point[i], 3), "ci95": ci(i)} for i, k in enumerate(("precision", "recall", "f1"))}}


def main(variants):
    split = {p["paper_id"]: p["split"] for p in json.loads((HERE / "sample_all.json").read_text(encoding="utf-8"))["papers"]}
    gold_all = json.loads((HERE / "gold_all.json").read_text(encoding="utf-8"))
    gold = {p: g for p, g in gold_all.items() if split[p] == "Test"}
    rest = json.loads((HERE / "gold_test_rest.json").read_text(encoding="utf-8"))
    sets = {**{p: "biored_all" for p in gold}, **{p: "biored_testrest" for p in rest}}
    gold.update(rest)
    report, all_counts, all_failed = {}, {}, {}
    for v in variants:
        counts, failed = {}, []
        for p, run_set in sets.items():
            out = RUNS / run_set / v / p / "outputs"
            if (out / "validation.json").exists() or (out / "source.json").exists():
                counts[p] = score_paper(p, gold[p], out)[0]
            else:
                failed.append(p)
                counts[p] = dict(EMPTY, fn=len(gold[p]["relations"]))
        valid = sorted(p for p in counts if p not in failed)
        all_counts[v], all_failed[v] = counts, failed
        rng = random.Random(20261028)
        report[v] = {"abstracts": len(sets), "failed_or_missing": failed,
                     "valid_only": with_ci(counts, valid, False, rng),
                     "deployment": with_ci(counts, sorted(counts), False, rng),
                     "false_positives_on_27_without_relation": sum(counts[p]["our_supported"] for p in rest if p in valid)}
        r = report[v]
        print(f"{v}: valid {len(valid)}/{len(sets)}; valid-only P/R/F1 {r['valid_only']['precision']}/{r['valid_only']['recall']}/"
              f"{r['valid_only']['f1']} F1 CI {r['valid_only']['ci95']['f1']}; deployment F1 {r['deployment']['f1']}; "
              f"supported pairs on the 27 without relation: {r['false_positives_on_27_without_relation']}")
    report["paired"] = {f"{a} - {b}": paired(a, b, all_counts, all_failed)
                        for i, a in enumerate(variants) for b in variants[i + 1:]}
    print(json.dumps(report["paired"], indent=1))
    (RUNS / "biored_testrest" / "test100.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1:] or ["v1_card_lean_medium", "v6_minimal_lean_medium", "pubtator3"])
