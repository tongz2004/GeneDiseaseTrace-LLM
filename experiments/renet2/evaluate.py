"""Score a single-shot run on the RENET2 full-text sample (gold from prepare.py).

Usage: python experiments/renet2/evaluate.py [variant] [run_set]      (defaults: v1_card_lean_medium renet2)

Matching and counting are BioRED's (experiments/biored/evaluate.py, imported unchanged). Gold relations are the
pairs RENET2 labels associated (expert-confirmed). Two views of our predictions: supported, and supported +
uncertain. Because RENET2's negatives were never checked by hand and its positives started from model
predictions, each of our supported pairs outside the gold set is also classed as: labelled non-associated (all three
tools agreed there was nothing), labelled ambiguous, or never labelled at all (a pair no tool proposed).
A status table shows what our method said about every labelled pair. Gold pairs whose gene or disease does not occur
in the (truncated) text are listed separately: no reader of this text can find them.
Writes runs/<run_set>/<variant>/renet2_eval.json and renet2_status_rows.csv.
"""
import collections
import csv
import importlib.util
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent / "single_shot" / "runs"
_spec = importlib.util.spec_from_file_location("biored_evaluate", HERE.parent / "biored" / "evaluate.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


def subset(gold, labels):
    return {p: {**g, "relations": [r for r in g["relations"] if r["type"] in labels]} for p, g in gold.items()}


def main(variant="v1_card_lean_medium", run_set="renet2"):
    gold = json.loads((HERE / "gold_renet2.json").read_text(encoding="utf-8"))
    base = RUNS / run_set / variant
    papers = sorted(p for p in gold if (base / p / "outputs" / "validation.json").exists()
                    or (base / p / "outputs" / "source.json").exists())  # source.json: PubTator3 baseline
    out = lambda p: base / p / "outputs"
    rng = random.Random(20261028)
    in_text = {p: {m["identifier"] for m in gold[p]["mentions"]} for p in gold}
    unreachable = [(p, r["gene"], r["disease"]) for p in papers for r in gold[p]["relations"]
                   if r["type"] == "associated" and not {r["gene"], r["disease"]} <= in_text[p]]
    result = {"variant": variant, "run_set": run_set, "papers_scored": len(papers),
              "papers_missing": sorted(set(gold) - set(papers)), "associated_pairs_not_in_text": len(unreachable)}
    positives = subset(gold, {"associated"})
    for view, positive in (("supported", ("supported",)), ("supported_or_uncertain", ("supported", "uncertain"))):
        counts = {p: B.score_paper(p, positives[p], out(p), positive)[0] for p in papers}
        result[view] = {**B.with_ci(counts, papers, False, rng), "counts": B.total([counts[p] for p in papers])}
        # Phenotype-inclusive view: gold pairs whose disease our rules keep as a symptom/phenotype set aside.
        result[f"{view} | scope_adjusted"] = B.with_ci(counts, papers, True, rng)
    # Our supported pairs outside the gold set, by what RENET2 says about them.
    extra = collections.Counter()
    status_rows = []
    table = collections.defaultdict(collections.Counter)
    for p in papers:
        _, pos_rows = B.score_paper(p, positives[p], out(p), ("supported",))
        _, any_rows = B.score_paper(p, gold[p], out(p), ("supported",))
        labelled_hits = {r["our_pair"]: r["gold_type"] for r in any_rows if r["outcome"] == "TP"}
        for r in pos_rows:
            if r["outcome"] == "FP":
                extra[labelled_hits.get(r["our_pair"], "never labelled") if r["detail"] == "gold has both concepts but no relation"
                      else "duplicate of a credited pair" if r["detail"].startswith("duplicate")
                      else "entity not a RENET2 concept"] += 1
        for r in any_rows:
            if r["outcome"] == "FP":
                continue
            status = "supported" if r["outcome"] == "TP" else (r["detail"].replace("our pair: ", "") if r["detail"].startswith("our pair")
                                                              else "entity not recorded")
            table[r["gold_type"]][status] += 1
            status_rows.append({**r, "our_status_for_pair": status})
    result["supported_not_in_gold"] = dict(extra)
    result["status_by_gold_label"] = {k: dict(v) for k, v in sorted(table.items())}
    (base / "renet2_eval.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["paper_id", "gold_type", "gold_gene", "gold_disease", "outcome", "detail", "our_status_for_pair",
              "our_gene", "our_disease", "our_relation", "our_source", "our_reason"]
    with (base / "renet2_status_rows.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(status_rows)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
