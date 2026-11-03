"""Blind, mixed check of the full-text pairs V6 (Claude) and BioREx report but BioTriplex links no relation for.

Usage: python experiments/biotriplex/precision_compare.py            (build the sheet)
       python experiments/biotriplex/precision_compare.py --score    (after the sheet is filled in)

Pools (precision_check.supported_concepts, one concept pair per supported pair whose ends are both BioTriplex
entities but which BioTriplex links no relation for), on the articles both runs have a valid result for:
- V6: runs/biotriplex/v6_minimal_lean_medium, minus the 20 pairs already judged in precision_check/;
- BioREx: runs/biotriplex/biorex_fulltext.
15 drawn from each (random.Random(20261105)); a pair drawn for both is shown once and counted for both. Controls:
5 BioTriplex relations and 5 co-mentioned pairs no run supports and BioTriplex does not relate, none of them
used in precision_check/. Items are shuffled and numbered Q01 ...; every item shows paragraphs found with
BioTriplex's, V6's and BioREx's mentions alike (precision_check.shown_spans), and nothing shows which system
reported it. Writes precision_compare/sheet.html (export biotriplex_precision_compare_answers.csv) and key.csv;
--score writes result.json with, per system, how many of its drawn pairs the annotator found stated.
Exploratory: one annotator.
"""
import argparse
import collections
import csv
import json
import random
from pathlib import Path

import precision_check as pc

HERE = Path(__file__).resolve().parent
OUT = HERE / "precision_compare"
RUNS = HERE.parent / "single_shot" / "runs" / "biotriplex"
SYSTEMS = {"v6": RUNS / "v6_minimal_lean_medium", "biorex": RUNS / "biorex_fulltext"}
SEED, PER_SYSTEM, CONTROLS = 20261105, 15, 5
EXPORT = "biotriplex_precision_compare_answers.csv"


def build():
    gold_all = json.loads((HERE / "gold_biotriplex.json").read_text(encoding="utf-8"))
    papers = json.loads((HERE / "sample_biotriplex.json").read_text(encoding="utf-8"))["evaluation"]
    used = {(r["paper_id"], r["gene"], r["disease"]) for r in csv.DictReader((pc.OUT / "key.csv").open(encoding="utf-8-sig"))}
    pools, positives, negatives = {s: [] for s in SYSTEMS}, [], []
    for paper in papers:
        if not all((run / paper / "outputs" / ("validation.json" if s == "v6" else "source.json")).exists()
                   for s, run in SYSTEMS.items()):
            continue
        gold = gold_all[paper]
        spans = pc.concept_spans(gold)
        related = {(r["gene"], r["disease"]) for r in gold["relations"]}
        covered = set()
        for s, run in SYSTEMS.items():
            cov, unannotated = pc.supported_concepts(paper, spans, related, run)
            covered |= cov
            pools[s] += [(paper, g, d) for g, d in sorted(unannotated) if (paper, g, d) not in used]
        positives += [(paper, g, d) for g, d in sorted(related) if (paper, g, d) not in used]
        blocks_of = lambda c: {b for b, _, _ in spans[c]}
        negatives += [(paper, g, d) for g in sorted(c for c in spans if c.startswith("G"))
                      for d in sorted(c for c in spans if c.startswith("D"))
                      if (g, d) not in related and (g, d) not in covered and blocks_of(g) & blocks_of(d)
                      and (paper, g, d) not in used]
    rng = random.Random(SEED)
    drawn = collections.defaultdict(set)
    for s in SYSTEMS:
        for x in rng.sample(pools[s], PER_SYSTEM):
            drawn[x].add(s)
    items = [(x, "+".join(sorted(s))) for x, s in drawn.items()]
    items += [(x, "positive") for x in rng.sample(positives, CONTROLS)] + [(x, "negative") for x in rng.sample(negatives, CONTROLS)]
    rng.shuffle(items)
    rows, key = [], []
    for n, ((paper, g, d), kind) in enumerate(items, 1):
        gold, spans = gold_all[paper], pc.concept_spans(gold_all[paper])
        article = pc.parse_article(HERE / "papers" / f"{paper}.xml")
        gs = pc.shown_spans(paper, spans, g, tuple(SYSTEMS.values()))
        ds = pc.shown_spans(paper, spans, d, tuple(SYSTEMS.values()))
        gb, db = {b for b, _, _ in gs}, {b for b, _, _ in ds}
        both = [b for b in article["blocks"] if b["block_id"] in gb & db][:pc.MAX_PARAGRAPHS]
        shown = both or ([b for b in article["blocks"] if b["block_id"] in gb][:3]
                         + [b for b in article["blocks"] if b["block_id"] in db][:3])
        body = "".join(f'<p style="color:#666;font-size:.85em;margin-bottom:0">{b["section"] or "Abstract"}</p>'
                       + pc.marked_text([b], gs, ds) for b in sorted(shown, key=lambda b: b["block_id"]))
        texts = gold["entity_texts"]
        item_id = f"Q{n:02d}"
        rows.append({"item_id": item_id, "gene": " / ".join(texts[g][:3]), "disease": " / ".join(texts[d][:3]),
                     "abstract_html": body})
        key.append({"item_id": item_id, "kind": kind, "paper_id": paper, "gene": g, "disease": d,
                    "paragraphs": "both" if both else "separate",
                    **{f"pool_{s}": len(pools[s]) for s in SYSTEMS}})
    OUT.mkdir(exist_ok=True)
    pc.write_form(OUT / "sheet.html", rows, scope="这些段落", store="biotriplex_precision_compare_v1", download=EXPORT)
    with (OUT / "key.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(key[0]))
        writer.writeheader()
        writer.writerows(key)
    print(f"{len(rows)} items -> {OUT}; pools {{{', '.join(f'{s}: {len(p)}' for s, p in pools.items())}}}; "
          f"{collections.Counter(k['kind'] for k in key)}")


def score(path):
    if not path:
        path = next(p for p in (OUT / EXPORT, Path.home() / "Downloads" / EXPORT) if p.exists())
    print(f"answers: {path}")
    key = {r["item_id"]: r for r in csv.DictReader((OUT / "key.csv").open(encoding="utf-8-sig"))}
    report = collections.defaultdict(lambda: {"items": 0, "states": collections.Counter()})
    for a in csv.DictReader(Path(path).open(encoding="utf-8-sig")):
        state = a[pc.ROUND1[1]].split(":")[0].strip().lower()
        if state not in pc.STATES:
            print(f"{a['item_id']}: unreadable answer {a[pc.ROUND1[1]]!r}")
            continue
        kind = key[a["item_id"]]["kind"]
        for group in (kind.split("+") if kind not in ("positive", "negative") else [kind]):
            report[group]["items"] += 1
            report[group]["states"][state] += 1
    for e in report.values():
        e["relation_or_hedged"] = e["states"]["relation"] + e["states"]["hedged relation"]
        e["states"] = dict(e["states"])
    (OUT / "result.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--score", nargs="?", const="", metavar="ANSWERS_CSV")
    args = parser.parse_args()
    score(args.score) if args.score is not None else build()
