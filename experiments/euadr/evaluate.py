"""Score a single-shot run on the EU-ADR gene–disease subset (gold from prepare.py).

Usage: python experiments/euadr/evaluate.py [variant] [run_set]      (defaults: v1_card_lean_medium euadr)

Matching and counting are BioRED's (experiments/biored/evaluate.py, imported unchanged): entities aligned by
overlapping text positions, one gene–disease concept pair per abstract, bootstrap 95% intervals over abstracts.
EU-ADR adds labels, so two scopes are scored:
- PA:    positive associations are the gold relations;
- PA+SA: speculative associations count as relations too.
Each scope is scored with our supported pairs as predictions, and with supported + uncertain.
Status table: for every gold pair (PA, SA, NA), what our method said about it (supported / uncertain /
explicit_negative / not_found, or an entity we did not record). This tests the uncertain and negative states,
which BioRED cannot. Writes runs/<run_set>/<variant>/euadr_eval.json and euadr_status_rows.csv.
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


def status_of(row):
    if row["outcome"] == "TP":
        return "supported"
    detail = row["detail"]
    return detail.replace("our pair: ", "") if detail.startswith("our pair") else "entity not recorded"


def main(variant="v1_card_lean_medium", run_set="euadr"):
    gold = json.loads((HERE / "gold_euadr.json").read_text(encoding="utf-8"))
    base = RUNS / run_set / variant
    papers = sorted(p for p in gold if (base / p / "outputs" / "validation.json").exists()
                    or (base / p / "outputs" / "source.json").exists())  # source.json: PubTator3 baseline
    out = lambda p: base / p / "outputs"
    rng = random.Random(20261028)
    result = {"variant": variant, "run_set": run_set, "papers_scored": len(papers),
              "papers_missing": sorted(set(gold) - set(papers))}
    for scope, labels in (("PA", {"PA"}), ("PA+SA", {"PA", "SA"})):
        g = subset(gold, labels)
        for view, positive in (("supported", ("supported",)), ("supported_or_uncertain", ("supported", "uncertain"))):
            counts = {p: B.score_paper(p, g[p], out(p), positive)[0] for p in papers}
            result[f"{scope} | {view}"] = {**B.with_ci(counts, papers, False, rng), "counts": B.total([counts[p] for p in papers])}
            # Phenotype-inclusive view: gold pairs whose disease our rules keep as a symptom/phenotype set aside.
            result[f"{scope} | {view} | scope_adjusted"] = B.with_ci(counts, papers, True, rng)
    # Status table over all gold pairs, whatever their label.
    rows = []
    for p in papers:
        _, paper_rows = B.score_paper(p, gold[p], out(p), ("supported",))
        rows += [r for r in paper_rows if r["outcome"] != "FP"]
    table = collections.defaultdict(collections.Counter)
    for r in rows:
        table[r["gold_type"]][status_of(r)] += 1
    result["status_by_gold_label"] = {k: dict(v) for k, v in sorted(table.items())}
    (base / "euadr_eval.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["paper_id", "gold_type", "gold_gene", "gold_disease", "outcome", "detail", "our_status", "our_relation", "our_reason"]
    with (base / "euadr_status_rows.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
