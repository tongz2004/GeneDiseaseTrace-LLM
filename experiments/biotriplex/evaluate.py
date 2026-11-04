"""Score a single-shot run on BioTriplex (gold from prepare.py).

Usage: python experiments/biotriplex/evaluate.py [variant] [run_set] [--abstract]
       defaults: v6_minimal_lean_medium biotriplex; --abstract scores the BTXA (abstract-only) articles

Matching and counting are BioRED's (experiments/biored/evaluate.py, imported unchanged, one-to-one). The gold set is
the same for BTX and BTXA articles, so an abstract-only run is charged with every relation in the article.
Reported for the evaluation set fixed in sample_biotriplex.json (missing or failed articles are listed, and a
deployment view counts them as empty answers):
- P/R/F1 with paper-bootstrap intervals, for supported and for supported + uncertain;
- recall split by where BioTriplex links the relation: in the abstract (at least one link there) or only outside it;
- source agreement on true positives: our stated source (this study / background / both) against the locations
  (result / background) of BioTriplex's links for that pair.
Writes runs/<run_set>/<variant>/biotriplex_eval[_abstract].json and the matching _rows.csv.
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
SOURCE = {"this study": {"result"}, "background": {"background"}, "both": {"result", "background"}}
EMPTY = {"tp_gold": 0, "fn": 0, "fn_symptom": 0, "our_supported": 0, "tp_ours": 0, "uncertain_or_negative": 0,
         "uncertain_or_negative_in_gold": 0, "tp_ours_many_to_one": 0, "tp_gold_many_to_one": 0}


def main(variant="v6_minimal_lean_medium", run_set="biotriplex", abstract=False):
    suffix = "2" if run_set == "biotriplex2" else ""  # second validation batch, see batch2_protocol.md
    gold = json.loads((HERE / f"gold_biotriplex{suffix}.json").read_text(encoding="utf-8"))
    sample = json.loads((HERE / f"sample_biotriplex{suffix}.json").read_text(encoding="utf-8"))
    papers = sample["evaluation_abstract" if abstract else "evaluation"]
    base = RUNS / run_set / variant
    out = lambda p: base / p / "outputs"
    done = [p for p in papers if (out(p) / "validation.json").exists() or (out(p) / "source.json").exists()]
    result = {"variant": variant, "run_set": run_set, "input": "abstract only" if abstract else "full text",
              "scorer": B.SCORER, "papers": len(papers), "papers_scored": len(done),
              "papers_missing": sorted(set(papers) - set(done)), "gold_relations": sum(len(gold[p]["relations"]) for p in papers)}
    rows_out, where, source = [], collections.Counter(), collections.Counter()
    for view, positive in (("supported", ("supported",)), ("supported_or_uncertain", ("supported", "uncertain"))):
        counts = {}
        for p in papers:
            if p not in done:
                counts[p] = dict(EMPTY, fn=len(gold[p]["relations"]))
                continue
            counts[p], rows = B.score_paper(p, gold[p], out(p), positive)
            if view != "supported":
                continue
            for r, row in zip(gold[p]["relations"], rows):
                assert row["outcome"] in ("TP", "FN")
                place = "only outside the abstract" if r["body_only"] else "in the abstract"
                where[(place, row["outcome"])] += 1
                row.update(gold_locations="+".join(r["locations"]), gold_place=place, gold_types="; ".join(r["relation_types"]))
                if row["outcome"] == "TP":
                    ours, theirs = SOURCE.get(row["our_source"]), set(r["locations"]) - {""}
                    source["no source given" if ours is None else "same" if ours == theirs
                           else "overlapping" if ours & theirs else "different"] += 1
            rows_out += rows
        rng = random.Random(20261028)
        result[view] = {"valid_only": B.with_ci(counts, done, False, rng), "deployment": B.with_ci(counts, papers, False, rng),
                        "counts": B.total([counts[p] for p in done])}
    missing = [p for p in papers if p not in done]
    for place in ("in the abstract", "only outside the abstract"):
        tp, fn = where[(place, "TP")], where[(place, "FN")]
        result[f"recall | relations {place}"] = {"gold": tp + fn, "found": tp,
                                                 "recall": round(tp / (tp + fn), 3) if tp + fn else None}
    result["source_agreement_on_true_positives"] = dict(source)
    result["note"] = ("valid_only scores articles with a valid result; deployment counts the missing ones "
                      f"({len(missing)}) as empty answers. Recall by place and source agreement use valid articles.")
    stem = "biotriplex_eval_abstract" if abstract else "biotriplex_eval"
    (base / f"{stem}.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["paper_id", "outcome", "detail", "gold_gene", "gold_disease", "gold_place", "gold_locations", "gold_types",
              "our_pair", "our_gene", "our_disease", "our_status", "our_relation", "our_source", "our_reason"]
    with (base / f"{stem}_rows.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--abstract"]
    main(*args[:2], abstract="--abstract" in sys.argv)
